import os
import re
import uuid
import json
import asyncio
from datetime import datetime, timedelta
from io import BytesIO

from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from werkzeug.utils import secure_filename
from sqlalchemy import text, inspect
from pypdf import PdfReader
from deep_translator import GoogleTranslator

# IMPORTAÇÕES LOCAIS
from models import db, Deck, Card, Story, Sentence, StudyLog, Reference, Project, Draft
import utils

app = Flask(__name__)

# --- CONFIGURAÇÃO DO BANCO DE DADOS ---
app.config['SQLALCHEMY_DATABASE_URI'] = os.environ.get('DATABASE_URL', 'sqlite:///scifluency.db')
if app.config['SQLALCHEMY_DATABASE_URI'] and app.config['SQLALCHEMY_DATABASE_URI'].startswith("postgres://"):
    app.config['SQLALCHEMY_DATABASE_URI'] = app.config['SQLALCHEMY_DATABASE_URI'].replace("postgres://", "postgresql://", 1)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# --- CONFIGURAÇÃO DE UPLOAD ---
UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

db.init_app(app)
APP_NAME = "SciFluency Ultimate"

# --- AUTO-REPARO DO BANCO DE DADOS ---
def fix_db():
    """Garante que colunas novas existam sem apagar dados antigos"""
    with app.app_context():
        try:
            db.create_all()
            inspector = inspect(db.engine)
            
            # Verifica tabela CARD (para fonética e SRS)
            if inspector.has_table('card'):
                cols = [c['name'] for c in inspector.get_columns('card')]
                with db.engine.connect() as conn:
                    if 'ease_factor' not in cols: 
                        conn.execute(text("ALTER TABLE card ADD COLUMN ease_factor FLOAT DEFAULT 2.5"))
                        conn.commit()
                    if 'context' not in cols: 
                        conn.execute(text("ALTER TABLE card ADD COLUMN context TEXT"))
                        conn.commit()
                    if 'ipa' not in cols:
                        conn.execute(text("ALTER TABLE card ADD COLUMN ipa VARCHAR(200)"))
                        conn.commit()
                        
        except Exception as e:
            print(f"Erro no auto-reparo: {e}")

fix_db()

# --- ROTA INICIAL (DASHBOARD) ---
@app.route('/')
def index():
    # Inicializa dados básicos se vazio
    if not Project.query.first():
        try: db.session.add(Project(id="thesis", title="Meu Projeto", target_journal="A Definir")); db.session.commit()
        except: db.session.rollback()
    
    if not Deck.query.get("my_vocab"):
        try: db.session.add(Deck(id="my_vocab", name="Vocabulário", icon="🎓")); db.session.commit()
        except: pass
    
    # Estatísticas
    stats = {
        "vocab": Card.query.count(),
        "refs": Reference.query.count(),
        "stories": Story.query.count()
    }
    
    return render_template('layout.html', mode='dashboard', stats=stats, app_name=APP_NAME)

# --- FERRAMENTAS RESTAURADAS (PUBMED, CHECKER, ETC) ---

@app.route('/search', methods=['GET', 'POST'])
def search():
    results = []
    if request.method == 'POST':
        query = request.form.get('query')
        results = utils.search_pubmed(query)
    return render_template('layout.html', mode='search', results=results, app_name=APP_NAME)

@app.route('/checker', methods=['GET', 'POST'])
def checker():
    original = ""
    corrected = ""
    if request.method == 'POST':
        original = request.form.get('text_input')
        corrected = utils.improve_english_text(original)
    return render_template('layout.html', mode='checker', original=original, corrected=corrected, app_name=APP_NAME)

@app.route('/miner', methods=['GET', 'POST'])
def miner():
    keywords = []
    if request.method == 'POST':
        text = request.form.get('text_input', '')
        f = request.files.get('file_input')
        if f and f.filename.endswith('.pdf'):
            try: text += " ".join([p.extract_text() for p in PdfReader(f).pages])
            except: pass
        keywords = utils.get_top_keywords(text)
    return render_template('layout.html', mode='miner', keywords=keywords, app_name=APP_NAME)

@app.route('/summarizer', methods=['GET', 'POST'])
def summarizer():
    if request.method == 'POST':
        text = request.form.get('text_input', '')
        f = request.files.get('file_input')
        if f: 
            try: text = " ".join([p.extract_text() for p in PdfReader(f).pages])
            except: pass
        
        # Cria lição de resumo
        sid = str(uuid.uuid4())[:8]
        db.session.add(Story(id=sid, title="Resumo: " + text[:30] + "..."))
        
        # Usa formatação inteligente e limpa HTML para salvar frases
        formatted = utils.format_abstract_smart(text)
        clean = re.sub(r'<[^>]*>', '', formatted)
        
        for s in re.split(r'(?<=[.!?])\s+', clean)[:50]:
            if len(s) > 5: db.session.add(Sentence(en=s, pt="...", story_id=sid))
            
        db.session.commit()
        return redirect(url_for('ler', id=sid))
        
    return render_template('layout.html', mode='summarizer', app_name=APP_NAME)

# --- ROTAS DE ESTUDO (KARAOKÊ E TEXTO) ---

@app.route('/novo', methods=['GET', 'POST'])
def novo():
    if request.method == 'POST':
        text = request.form.get('text_input', '')
        f = request.files.get('file_input')
        
        # Extrai PDF
        if f: 
            try: text = " ".join([p.extract_text() for p in PdfReader(f).pages])
            except: pass
        
        if text:
            sid = str(uuid.uuid4())[:8]
            db.session.add(Story(id=sid, title=text[:40]+"..."))
            
            # Traduz e salva frases
            tr = GoogleTranslator(source='en', target='pt')
            # Limpa espaços e quebra
            clean_text = re.sub(r'\s+', ' ', text)
            for s in re.split(r'(?<=[.!?])\s+', clean_text)[:100]:
                if len(s) > 3:
                    try: pt = tr.translate(s)
                    except: pt="..."
                    db.session.add(Sentence(en=s, pt=pt, story_id=sid))
            
            db.session.commit()
            return redirect(url_for('ler', id=sid))
            
    return render_template('layout.html', mode='new', app_name=APP_NAME)

@app.route('/ler/<id>')
def ler(id):
    story = Story.query.get(id)
    return render_template('layout.html', mode='read', story=story, app_name=APP_NAME)

@app.route('/pronunciation')
def pronunciation():
    # Pega palavras aleatórias para treinar
    cards = Card.query.order_by(db.func.random()).limit(10).all()
    words = [{"w": c.front, "ipa": c.ipa} for c in cards]
    if not words: words = [{"w":"Science", "ipa":"/ˈsaɪ.əns/"}]
    return render_template('layout.html', mode='pronunciation', words_json=json.dumps(words), app_name=APP_NAME)

# --- ROTAS DE PESQUISA (BIBLIOTECA PDF) ---

@app.route('/library', methods=['GET', 'POST'])
def library():
    if request.method == 'POST':
        f = request.files.get('pdf_file')
        if f:
            fn = secure_filename(f.filename)
            f.save(os.path.join(app.config['UPLOAD_FOLDER'], fn))
            db.session.add(Reference(title=fn, status='to_read', pdf_filename=fn, project_id="thesis"))
            db.session.commit()
    return render_template('layout.html', mode='library', references=Reference.query.all(), app_name=APP_NAME)

@app.route('/read_ref/<int:id>')
def read_ref(id):
    # Converte PDF da biblioteca acadêmica em Lição de Estudo
    ref = Reference.query.get(id)
    try:
        path = os.path.join(app.config['UPLOAD_FOLDER'], ref.pdf_filename)
        txt = " ".join([p.extract_text() for p in PdfReader(path).pages])
        
        sid = str(uuid.uuid4())[:8]
        db.session.add(Story(id=sid, title="PDF: "+ref.title))
        
        for s in re.split(r'(?<=[.!?])\s+', re.sub(r'\s+', ' ', txt))[:150]:
            if len(s) > 5: db.session.add(Sentence(en=s, pt="...", story_id=sid))
            
        db.session.commit()
        return redirect(url_for('ler', id=sid))
    except: return "Erro ao ler PDF"

@app.route('/writer')
def writer():
    return render_template('layout.html', mode='writer', connectors=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

# --- ROTAS AUXILIARES E VOCABULÁRIO (COM FONÉTICA) ---

@app.route('/vocab_list')
def vocab_list():
    return render_template('layout.html', mode='vocab_list', cards=Card.query.all(), app_name=APP_NAME)

@app.route('/traduzir_palavra')
def traduzir():
    w = request.args.get('w', '')
    # Traduz
    try: t = GoogleTranslator(source='en', target='pt').translate(w)
    except: t = "..."
    
    # Gera IPA
    ipa_text = utils.get_phonetic(w)
    
    return jsonify({"t": t, "ipa": ipa_text})

@app.route('/adicionar_vocab', methods=['POST'])
def add_vocab():
    term = request.form.get('term')
    
    # Gera IPA antes de salvar
    ipa_text = utils.get_phonetic(term)
    
    db.session.add(Card(
        front=term, 
        back="...", 
        ipa=ipa_text,  # Salva no banco
        deck_id="my_vocab", 
        next_review=datetime.now().strftime('%Y-%m-%d')
    ))
    
    db.session.commit()
    return jsonify({'status':'ok'})

@app.route('/tts', methods=['POST'])
def tts():
    audio = utils.get_audio_sync(request.form.get('text'), request.form.get('accent'))
    return send_file(audio, mimetype='audio/mpeg') if audio else ("Error", 500)

@app.route('/delete_story/<id>')
def delete_story(id):
    db.session.delete(Story.query.get(id)); db.session.commit()
    return redirect(url_for('index'))

@app.route('/delete_ref/<int:id>')
def delete_ref(id):
    db.session.delete(Reference.query.get(id)); db.session.commit()
    return redirect(url_for('library'))

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
