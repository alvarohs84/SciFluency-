import os
import re
import json
import uuid
from datetime import datetime, timedelta
from io import BytesIO

from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from werkzeug.utils import secure_filename
from sqlalchemy import text, inspect
from pypdf import PdfReader
from deep_translator import GoogleTranslator
from gtts import gTTS
import edge_tts
import asyncio

# IMPORTAÇÕES LOCAIS
from models import db, Deck, Card, Story, Sentence, StudyLog, Reference, Project, Note, Draft
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
APP_NAME = "SciFluency Ultimate"

# --- AUTO-REPARO DO BANCO ---
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

# --- ROTAS PRINCIPAIS ---

@app.route('/')
def index():
    if not Project.query.first():
        try:
            db.session.add(Project(id="thesis", title="Meu Projeto", target_journal="Science"))
            if not Deck.query.get("my_vocab"): db.session.add(Deck(id="my_vocab", name="Vocabulário", icon="🎓"))
            db.session.commit()
        except: db.session.rollback()
    
    # Stats combinados
    stats = {
        "vocab": Card.query.count(),
        "refs": Reference.query.count(),
        "heatmap": []
    }
    for i in range(6, -1, -1):
        d = (datetime.now() - timedelta(days=i)).strftime('%Y-%m-%d')
        l = StudyLog.query.filter_by(date=d).first()
        stats['heatmap'].append({"date": d, "color": "#00b894" if l else "#ecf0f1"})

    return render_template('layout.html', mode='dashboard', project=Project.query.first(), stats=stats, app_name=APP_NAME)

# --- MÓDULO 1: LEITURA CLÁSSICA (TEXTO) ---
@app.route('/novo_texto')
def novo_texto():
    return render_template('layout.html', mode='new_text', app_name=APP_NAME)

@app.route('/processar', methods=['POST'])
def processar():
    text_content = request.form.get('texto_full', '')
    if text_content.strip():
        sid = str(uuid.uuid4())[:8]
        # Limpa e divide em frases
        clean_text = re.sub(r'\s+', ' ', text_content)
        frases = re.split(r'(?<=[.!?])\s+(?=[A-Z])', clean_text)
        
        # Cria a Story
        title = "Leitura: " + (frases[0][:30] + "..." if frases else "Texto Manual")
        db.session.add(Story(id=sid, title=title))
        
        # Traduz e salva frases
        tr = GoogleTranslator(source='en', target='pt')
        for f in frases[:50]: # Limite de 50 frases por performance
            if len(f) < 5: continue
            try: pt = tr.translate(f)
            except: pt = "..."
            db.session.add(Sentence(en=f, pt=pt, story_id=sid))
        
        db.session.commit()
        return redirect(url_for('ler', id=sid))
    return redirect(url_for('index'))

@app.route('/ler/<id>')
def ler(id):
    story = Story.query.get(id)
    return render_template('layout.html', mode='read_classic', story=story, app_name=APP_NAME)

# --- MÓDULO 2: PESQUISA (PDFs) ---
@app.route('/library', methods=['GET', 'POST'])
def library():
    if request.method == 'POST':
        f = request.files.get('pdf_file')
        if f:
            filename = secure_filename(f.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            f.save(filepath)
            try: abstract = PdfReader(filepath).pages[0].extract_text()[:400]
            except: abstract = "..."
            db.session.add(Reference(title=filename, status='to_read', pdf_filename=filename, abstract=abstract, project_id="thesis"))
            db.session.commit()
            return redirect(url_for('library'))
    return render_template('layout.html', mode='library', references=Reference.query.all(), app_name=APP_NAME)

@app.route('/read_ref/<int:id>')
def read_ref(id):
    ref = Reference.query.get(id)
    try:
        path = os.path.join(app.config['UPLOAD_FOLDER'], ref.pdf_filename)
        text_content = " ".join([p.extract_text() or "" for p in PdfReader(path).pages])
        formatted = utils.format_abstract_smart(text_content)
        return render_template('layout.html', mode='read_pdf', title=ref.title, content=formatted, app_name=APP_NAME)
    except: return "Erro ao ler PDF"

# --- MÓDULO 3: FONÉTICA & SPEAKING ---
@app.route('/pronunciation')
def pronunciation():
    # Pega palavras do deck para treinar
    cards = Card.query.order_by(db.func.random()).limit(10).all()
    words_data = [{"w": c.front, "ipa": c.ipa} for c in cards]
    if not words_data: words_data = [{"w": "Science", "ipa": "/ˈsaɪ.əns/"}]
    return render_template('layout.html', mode='pronunciation', words_json=json.dumps(words_data), app_name=APP_NAME)

# --- ROTAS UTILITÁRIAS ---
@app.route('/vocab_list')
def vocab_list(): return render_template('layout.html', mode='vocab_list', cards=Card.query.all(), app_name=APP_NAME)

@app.route('/writer')
def writer(): return render_template('layout.html', mode='writer', connectors=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

@app.route('/tts', methods=['POST'])
def tts():
    audio = utils.get_audio_sync(request.form.get('text'), request.form.get('accent'))
    return send_file(audio, mimetype='audio/mpeg') if audio else ("Error", 500)

@app.route('/traduzir_palavra')
def traduzir(): return jsonify({"t": GoogleTranslator(source='en', target='pt').translate(request.args.get('w'))})

@app.route('/adicionar_vocab', methods=['POST'])
def add_vocab():
    db.session.add(Card(front=request.form.get('term'), back="...", deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d')))
    db.session.commit()
    return jsonify({"status": "ok"})

@app.route('/delete_ref/<int:id>')
def delete_ref(id):
    r = Reference.query.get(id)
    db.session.delete(r); db.session.commit()
    return redirect(url_for('library'))

if __name__ == '__main__': app.run(host='0.0.0.0', port=5000)
