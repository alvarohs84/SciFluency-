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
from gtts import gTTS
import edge_tts

from models import db, Deck, Card, Story, Sentence, StudyLog, Reference, Project, Draft
import utils

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = os.environ.get('DATABASE_URL', 'sqlite:///scifluency.db')
if app.config['SQLALCHEMY_DATABASE_URI'].startswith("postgres://"):
    app.config['SQLALCHEMY_DATABASE_URI'] = app.config['SQLALCHEMY_DATABASE_URI'].replace("postgres://", "postgresql://", 1)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

db.init_app(app)
APP_NAME = "SciFluency Pro"

# --- REPARO AUTOMÁTICO ---
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
        except: pass
fix_db()

# --- ROTAS DE APRENDIZADO (FUNCIONALIDADES ANTIGAS RESTAURADAS) ---

@app.route('/')
def index():
    if not Project.query.first():
        try:
            db.session.add(Project(id="thesis", title="Meu Trabalho", target_journal="A Definir"))
            if not Deck.query.get("my_vocab"): db.session.add(Deck(id="my_vocab", name="Vocabulário", icon="🎓"))
            db.session.commit()
        except: db.session.rollback()
    
    # Estatísticas completas
    stats = {
        "vocab": Card.query.count(),
        "stories": Story.query.count(),
        "refs": Reference.query.count(),
        "heatmap": []
    }
    for i in range(6, -1, -1):
        d = (datetime.now() - timedelta(days=i)).strftime('%Y-%m-%d')
        l = StudyLog.query.filter_by(date=d).first()
        stats['heatmap'].append({"date": d, "color": "#00b894" if l else "#ecf0f1"})

    return render_template('layout.html', mode='dashboard', stats=stats, project=Project.query.first(), app_name=APP_NAME)

@app.route('/processar', methods=['POST'])
def processar():
    """O CORAÇÃO DO SISTEMA: Pega texto/PDF e transforma em Lição Interativa"""
    text_content = request.form.get('texto_full', '')
    uploaded_files = request.files.getlist("arquivo_upload")
    
    # Extrai de PDF se houver
    for f in uploaded_files:
        if f and f.filename.endswith('.pdf'):
            try: text_content += " ".join([p.extract_text() for p in PdfReader(f).pages])
            except: pass

    if text_content.strip():
        sid = str(uuid.uuid4())[:8]
        clean_text = re.sub(r'\s+', ' ', text_content)
        # Quebra inteligente de frases
        frases = re.split(r'(?<=[.!?])\s+(?=[A-Z])', clean_text)
        
        # Salva como 'Story' para poder usar o Karaokê
        title = "Estudo: " + (frases[0][:40] + "..." if frases else "Novo Texto")
        db.session.add(Story(id=sid, title=title))
        
        tr = GoogleTranslator(source='en', target='pt')
        for f in frases[:100]: # Limite para não travar
            if len(f) < 5: continue
            try: pt = tr.translate(f)
            except: pt = "..."
            db.session.add(Sentence(en=f, pt=pt, story_id=sid))
        
        db.session.commit()
        return redirect(url_for('ler', id=sid))
        
    return redirect(url_for('novo'))

@app.route('/novo')
def novo(): return render_template('layout.html', mode='new', app_name=APP_NAME)

@app.route('/ler/<id>')
def ler(id):
    story = Story.query.get(id)
    return render_template('layout.html', mode='read_interactive', story=story, app_name=APP_NAME)

@app.route('/pronunciation')
def pronunciation():
    cards = Card.query.order_by(db.func.random()).limit(10).all()
    words = [{"w": c.front, "ipa": c.ipa} for c in cards]
    if not words: words = [{"w": "Science", "ipa": "/ˈsaɪ.əns/"}]
    return render_template('layout.html', mode='pronunciation', words_json=json.dumps(words), app_name=APP_NAME)

@app.route('/vocab_list')
def vocab_list(): return render_template('layout.html', mode='vocab_list', cards=Card.query.all(), app_name=APP_NAME)

@app.route('/library_list')
def library_list():
    # Lista de textos processados (Stories)
    return render_template('layout.html', mode='library_list', stories=Story.query.all(), app_name=APP_NAME)

# --- ROTAS DE PESQUISA (FUNCIONALIDADES NOVAS ADICIONADAS) ---

@app.route('/ref_library', methods=['GET', 'POST'])
def ref_library():
    if request.method == 'POST':
        f = request.files.get('pdf_file')
        if f:
            filename = secure_filename(f.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            f.save(filepath)
            db.session.add(Reference(title=filename, status='to_read', pdf_filename=filename, project_id="thesis"))
            db.session.commit()
            return redirect(url_for('ref_library'))
    return render_template('layout.html', mode='ref_library', references=Reference.query.all(), app_name=APP_NAME)

@app.route('/read_pdf_ref/<int:id>')
def read_pdf_ref(id):
    # Lê PDF da biblioteca de pesquisa e joga no Leitor Interativo
    ref = Reference.query.get(id)
    try:
        path = os.path.join(app.config['UPLOAD_FOLDER'], ref.pdf_filename)
        text_content = " ".join([p.extract_text() or "" for p in PdfReader(path).pages])
        
        # Cria uma Story temporária ou permanente para leitura
        sid = str(uuid.uuid4())[:8]
        db.session.add(Story(id=sid, title="PDF: " + ref.title))
        frases = re.split(r'(?<=[.!?])\s+', text_content)
        for f in frases[:100]:
             if len(f)>5: db.session.add(Sentence(en=f, pt="...", story_id=sid))
        db.session.commit()
        return redirect(url_for('ler', id=sid))
    except: return "Erro ao ler PDF"

@app.route('/writer')
def writer(): return render_template('layout.html', mode='writer', connectors=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

# --- UTILITÁRIOS (TTS, Tradução, Save) ---
@app.route('/tts', methods=['POST'])
def tts():
    audio = utils.get_audio_sync(request.form.get('text'), request.form.get('accent'))
    return send_file(audio, mimetype='audio/mpeg') if audio else ("Error", 500)

@app.route('/traduzir_palavra')
def traduzir(): return jsonify({"t": GoogleTranslator(source='en', target='pt').translate(request.args.get('w'))})

@app.route('/adicionar_vocab', methods=['POST'])
def add_vocab():
    db.session.add(Card(front=request.form.get('term'), back="...", deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d')))
    
    # Log de atividade
    today = datetime.now().strftime('%Y-%m-%d')
    log = StudyLog.query.filter_by(date=today).first()
    if not log: log = StudyLog(date=today, count=0); db.session.add(log)
    log.count += 1
    
    db.session.commit()
    return jsonify({"status": "ok"})

@app.route('/delete_story/<id>')
def delete_story(id):
    s = Story.query.get(id)
    db.session.delete(s); db.session.commit()
    return redirect(url_for('library_list'))

if __name__ == '__main__': app.run(host='0.0.0.0', port=5000)
