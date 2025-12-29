import os
import uuid
import json
from datetime import datetime, timedelta
from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from werkzeug.utils import secure_filename
from sqlalchemy import text, inspect
from pypdf import PdfReader
from deep_translator import GoogleTranslator

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
        except: pass
fix_db()

@app.route('/')
def index():
    if not Project.query.first():
        try: db.session.add(Project(id="thesis", title="Meu Projeto", target_journal="Nature")); db.session.commit()
        except: db.session.rollback()
    if not Deck.query.get("my_vocab"):
        try: db.session.add(Deck(id="my_vocab", name="Vocabulário", icon="🎓")); db.session.commit()
        except: pass
    
    stats = {"vocab": Card.query.count(), "refs": Reference.query.count(), "stories": Story.query.count()}
    return render_template('layout.html', mode='dashboard', stats=stats, app_name=APP_NAME)

# --- ROTAS DE FERRAMENTAS RESTAURADAS ---

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
    # Versão simplificada: processa e joga pro leitor
    if request.method == 'POST':
        text = request.form.get('text_input', '')
        f = request.files.get('file_input')
        if f: 
            try: text = " ".join([p.extract_text() for p in PdfReader(f).pages])
            except: pass
        
        # Cria uma história temporária para leitura
        sid = str(uuid.uuid4())[:8]
        db.session.add(Story(id=sid, title="Summary: " + text[:30]))
        formatted = utils.format_abstract_smart(text) # Reusa a função de formatação
        # Precisamos salvar frases para o Karaoke funcionar
        clean = re.sub(r'<[^>]*>', '', formatted)
        for s in re.split(r'(?<=[.!?])\s+', clean)[:50]:
            if len(s)>5: db.session.add(Sentence(en=s, pt="...", story_id=sid))
        db.session.commit()
        return redirect(url_for('ler', id=sid))
    return render_template('layout.html', mode='summarizer', app_name=APP_NAME)

# --- ROTAS DE ESTUDO & PESQUISA (MANTIDAS) ---

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
def ler(id):
    return render_template('layout.html', mode='read', story=Story.query.get(id), app_name=APP_NAME)

@app.route('/pronunciation')
def pronunciation():
    words = [{"w": c.front, "ipa": c.ipa} for c in Card.query.limit(10).all()]
    if not words: words = [{"w":"Science", "ipa":"/ˈsaɪ.əns/"}]
    return render_template('layout.html', mode='pronunciation', words_json=json.dumps(words), app_name=APP_NAME)

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
    # Converte PDF da biblioteca em Lição de Karaokê
    ref = Reference.query.get(id)
    try:
        txt = " ".join([p.extract_text() for p in PdfReader(os.path.join(app.config['UPLOAD_FOLDER'], ref.pdf_filename)).pages])
        sid = str(uuid.uuid4())[:8]
        db.session.add(Story(id=sid, title="PDF: "+ref.title))
        for s in re.split(r'(?<=[.!?])\s+', re.sub(r'\s+', ' ', txt))[:150]:
            if len(s)>5: db.session.add(Sentence(en=s, pt="...", story_id=sid))
        db.session.commit()
        return redirect(url_for('ler', id=sid))
    except: return "Erro no PDF"

@app.route('/writer')
def writer(): return render_template('layout.html', mode='writer', connectors=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

@app.route('/vocab_list')
def vocab_list(): return render_template('layout.html', mode='vocab', cards=Card.query.all(), app_name=APP_NAME)

# Utilitários
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
    return jsonify({'status':'ok'})

@app.route('/delete_story/<id>')
def delete_story(id):
    db.session.delete(Story.query.get(id)); db.session.commit()
    return redirect(url_for('index'))

@app.route('/delete_ref/<int:id>')
def delete_ref(id):
    db.session.delete(Reference.query.get(id)); db.session.commit()
    return redirect(url_for('library'))

if __name__ == '__main__': app.run(host='0.0.0.0', port=5000)
