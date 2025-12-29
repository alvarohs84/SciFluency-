import os
import re
import uuid
import json
import asyncio
from datetime import datetime
from io import BytesIO

from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from werkzeug.utils import secure_filename
from sqlalchemy import text, inspect
from pypdf import PdfReader
from deep_translator import GoogleTranslator

from models import db, Deck, Card, Story, Sentence, StudyLog, Reference, Project
import utils

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = os.environ.get('DATABASE_URL', 'sqlite:///scifluency.db')
if app.config['SQLALCHEMY_DATABASE_URI'] and app.config['SQLALCHEMY_DATABASE_URI'].startswith("postgres://"):
    app.config['SQLALCHEMY_DATABASE_URI'] = app.config['SQLALCHEMY_DATABASE_URI'].replace("postgres://", "postgresql://", 1)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

db.init_app(app)
APP_NAME = "SciFluency Ultimate"

def fix_db():
    with app.app_context():
        try:
            db.create_all()
            inspector = inspect(db.engine)
            if inspector.has_table('card'):
                cols = [c['name'] for c in inspector.get_columns('card')]
                with db.engine.connect() as conn:
                    if 'ease_factor' not in cols: conn.execute(text("ALTER TABLE card ADD COLUMN ease_factor FLOAT DEFAULT 2.5")); conn.commit()
                    if 'context' not in cols: conn.execute(text("ALTER TABLE card ADD COLUMN context TEXT")); conn.commit()
                    if 'ipa' not in cols: conn.execute(text("ALTER TABLE card ADD COLUMN ipa VARCHAR(200)")); conn.commit()
        except: pass
fix_db()

@app.route('/')
def index():
    if not Project.query.first():
        try: db.session.add(Project(id="thesis", title="Meu Projeto", target_journal="Definir")); db.session.commit()
        except: db.session.rollback()
    if not Deck.query.get("my_vocab"):
        try: db.session.add(Deck(id="my_vocab", name="Vocabulário", icon="🎓")); db.session.commit()
        except: pass
    
    stats = {"vocab": Card.query.count(), "refs": Reference.query.count(), "stories": Story.query.count()}
    return render_template('layout.html', mode='dashboard', stats=stats, app_name=APP_NAME)

# --- UPLOAD INTELIGENTE (PDF, RIS, NBIB) ---
@app.route('/library', methods=['GET', 'POST'])
def library():
    if request.method == 'POST':
        # Pega LISTA de arquivos
        files = request.files.getlist('ref_files')
        
        for f in files:
            if not f or not f.filename: continue
            filename = secure_filename(f.filename)
            ext = filename.rsplit('.', 1)[1].lower() if '.' in filename else ''
            
            # CASO 1: É PDF (Salva arquivo, extrai texto se der)
            if ext == 'pdf':
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                f.save(filepath)
                
                # Tenta extrair título do PDF (Metadados básicos)
                title = filename
                try:
                    reader = PdfReader(filepath)
                    if reader.metadata and reader.metadata.title:
                        title = reader.metadata.title
                except: pass
                
                db.session.add(Reference(title=title, authors="Ver PDF", year="2024", status='to_read', pdf_filename=filename, project_id="thesis"))

            # CASO 2: É METADADOS (RIS, NBIB) - Lê conteúdo e cria referência
            elif ext in ['ris', 'nbib', 'txt']:
                content = f.read().decode('utf-8', errors='ignore')
                parsed_refs = utils.parse_bib_file(content)
                
                for r in parsed_refs:
                    # Cria referência sem PDF atrelado (apenas dados)
                    db.session.add(Reference(
                        title=r.get('title', 'Sem Título'),
                        authors=r.get('authors', 'Desconhecido'),
                        year=r.get('year', 's.d.'),
                        abstract=r.get('abstract', ''),
                        status='to_read',
                        pdf_filename="", # Sem PDF físico
                        project_id="thesis"
                    ))
        
        db.session.commit()
        return redirect(url_for('library'))
        
    return render_template('layout.html', mode='library', references=Reference.query.order_by(Reference.id.desc()).all(), app_name=APP_NAME)

@app.route('/read_ref/<int:id>')
def read_ref(id):
    ref = Reference.query.get(id)
    content = ""
    # Se tem PDF, lê o PDF
    if ref.pdf_filename:
        try:
            path = os.path.join(app.config['UPLOAD_FOLDER'], ref.pdf_filename)
            content = " ".join([p.extract_text() for p in PdfReader(path).pages])
        except: content = "Erro ao ler PDF."
    # Se não tem PDF (veio do RIS/NBIB), usa o Abstract
    else:
        content = ref.abstract if ref.abstract else "Resumo não disponível."
        
    sid = str(uuid.uuid4())[:8]
    db.session.add(Story(id=sid, title="Ref: "+ref.title))
    
    # Prepara Karaokê
    clean = re.sub(r'\s+', ' ', content)
    for s in re.split(r'(?<=[.!?])\s+', clean)[:150]:
        if len(s) > 5: db.session.add(Sentence(en=s, pt="...", story_id=sid))
        
    db.session.commit()
    return redirect(url_for('ler', id=sid))

# --- OUTRAS ROTAS (MANTIDAS) ---

@app.route('/novo', methods=['GET', 'POST'])
def novo():
    if request.method == 'POST':
        text = request.form.get('text_input', '')
        f = request.files.get('file_input')
        if f: 
            try: text = " ".join([p.extract_text() for p in PdfReader(f).pages])
            except: pass
        if text:
            sid = str(uuid.uuid4())[:8]
            db.session.add(Story(id=sid, title=text[:40]+"..."))
            tr = GoogleTranslator(source='en', target='pt')
            for s in re.split(r'(?<=[.!?])\s+', re.sub(r'\s+', ' ', text))[:100]:
                if len(s)>3:
                    try: pt = tr.translate(s)
                    except: pt="..."
                    db.session.add(Sentence(en=s, pt=pt, story_id=sid))
            db.session.commit()
            return redirect(url_for('ler', id=sid))
    return render_template('layout.html', mode='new', app_name=APP_NAME)

@app.route('/ler/<id>')
def ler(id): return render_template('layout.html', mode='read', story=Story.query.get(id), app_name=APP_NAME)

@app.route('/pronunciation')
def pronunciation():
    cards = Card.query.order_by(db.func.random()).limit(10).all()
    words = [{"w": c.front, "ipa": c.ipa} for c in cards]
    if not words: words = [{"w":"Science", "ipa":"/ˈsaɪ.əns/"}]
    return render_template('layout.html', mode='pronunciation', words_json=json.dumps(words), app_name=APP_NAME)

@app.route('/search', methods=['GET', 'POST'])
def search():
    results = []
    if request.method == 'POST':
        results = utils.search_pubmed(request.form.get('query'))
    return render_template('layout.html', mode='search', results=results, app_name=APP_NAME)

@app.route('/checker', methods=['GET', 'POST'])
def checker():
    orig, corr = "", ""
    if request.method == 'POST':
        orig = request.form.get('text_input')
        corr = utils.improve_english_text(orig)
    return render_template('layout.html', mode='checker', original=orig, corrected=corr, app_name=APP_NAME)

@app.route('/miner', methods=['GET', 'POST'])
def miner():
    k = []
    if request.method == 'POST':
        txt = request.form.get('text_input', '')
        f = request.files.get('file_input')
        if f: 
            try: txt += " ".join([p.extract_text() for p in PdfReader(f).pages])
            except: pass
        k = utils.get_top_keywords(txt)
    return render_template('layout.html', mode='miner', keywords=k, app_name=APP_NAME)

@app.route('/summarizer', methods=['GET', 'POST'])
def summarizer():
    if request.method == 'POST':
        txt = request.form.get('text_input', '')
        f = request.files.get('file_input')
        if f: 
            try: txt = " ".join([p.extract_text() for p in PdfReader(f).pages])
            except: pass
        sid = str(uuid.uuid4())[:8]
        db.session.add(Story(id=sid, title="Resumo: "+txt[:20]))
        clean = re.sub(r'<[^>]*>', '', utils.format_abstract_smart(txt))
        for s in re.split(r'(?<=[.!?])\s+', clean)[:50]:
            if len(s)>5: db.session.add(Sentence(en=s, pt="...", story_id=sid))
        db.session.commit()
        return redirect(url_for('ler', id=sid))
    return render_template('layout.html', mode='summarizer', app_name=APP_NAME)

@app.route('/writer')
def writer(): return render_template('layout.html', mode='writer', connectors=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

@app.route('/get_citation/<int:id>')
def get_citation(id):
    r = Reference.query.get(id)
    return jsonify(utils.generate_citation_formats(r)) if r else jsonify({"error": "404"})

@app.route('/vocab_list')
def vocab_list(): return render_template('layout.html', mode='vocab_list', cards=Card.query.all(), app_name=APP_NAME)

@app.route('/traduzir_palavra')
def traduzir():
    w = request.args.get('w','')
    try: t = GoogleTranslator(source='en', target='pt').translate(w)
    except: t="..."
    return jsonify({"t": t, "ipa": utils.get_phonetic(w)})

@app.route('/adicionar_vocab', methods=['POST'])
def add_vocab():
    t = request.form.get('term')
    db.session.add(Card(front=t, back="...", ipa=utils.get_phonetic(t), deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d')))
    db.session.commit()
    return jsonify({'status':'ok'})

@app.route('/tts', methods=['POST'])
def tts():
    a = utils.get_audio_sync(request.form.get('text'), request.form.get('accent'))
    return send_file(a, mimetype='audio/mpeg') if a else ("Erro",500)

@app.route('/delete_story/<id>')
def delete_story(id): db.session.delete(Story.query.get(id)); db.session.commit(); return redirect(url_for('index'))
@app.route('/delete_ref/<int:id>')
def delete_ref(id): db.session.delete(Reference.query.get(id)); db.session.commit(); return redirect(url_for('library'))

if __name__ == '__main__': app.run(host='0.0.0.0', port=5000)
