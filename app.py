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
from Bio import Entrez

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

# --- UPLOAD MULTIPLO & RIS ---
@app.route('/library', methods=['GET', 'POST'])
def library():
    if request.method == 'POST':
        files = request.files.getlist('ref_files')
        for f in files:
            if not f or not f.filename: continue
            filename = secure_filename(f.filename)
            ext = filename.rsplit('.', 1)[1].lower() if '.' in filename else ''
            
            if ext == 'pdf':
                f.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                title = filename
                try:
                    r = PdfReader(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                    if r.metadata and r.metadata.title: title = r.metadata.title
                except: pass
                db.session.add(Reference(title=title, authors="Ver PDF", year="2024", status='to_read', pdf_filename=filename, project_id="thesis"))

            elif ext in ['ris', 'nbib', 'txt']:
                parsed = utils.parse_bib_file(f.read().decode('utf-8', errors='ignore'))
                for r in parsed:
                    db.session.add(Reference(title=r.get('title','Sem Título'), authors=r.get('authors','Desconhecido'), year=r.get('year','s.d.'), abstract=r.get('abstract',''), status='to_read', pdf_filename="", project_id="thesis"))
        
        db.session.commit()
        return redirect(url_for('library'))
    return render_template('layout.html', mode='library', references=Reference.query.order_by(Reference.id.desc()).all(), app_name=APP_NAME)

# --- PUBMED IMPORT ---
@app.route('/import_pubmed/<pmid>')
def import_pubmed(pmid):
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
        raw = handle.read()
        title = re.search(r'TI  - (.*)', raw)
        title = title.group(1) if title else f"PubMed {pmid}"
        abstract = re.search(r'AB  - (.*)', raw, re.DOTALL)
        ab_clean = abstract.group(1)[:600]+"..." if abstract else "..."
        
        db.session.add(Reference(title=title, authors="Via PubMed", year=datetime.now().strftime('%Y'), status='to_read', pdf_filename="", abstract=ab_clean, project_id="thesis"))
        db.session.commit()
        return redirect(url_for('library'))
    except: return "Erro Import"

@app.route('/read_ref/<int:id>')
def read_ref(id):
    ref = Reference.query.get(id)
    content = ""
    if ref.pdf_filename:
        try: content = " ".join([p.extract_text() for p in PdfReader(os.path.join(app.config['UPLOAD_FOLDER'], ref.pdf_filename)).pages])
        except: content = "Erro PDF"
    else: content = ref.abstract if ref.abstract else "..."
    sid = str(uuid.uuid4())[:8]
    db.session.add(Story(id=sid, title="Ref: "+ref.title[:20]))
    clean = re.sub(r'\s+', ' ', content)
    for s in re.split(r'(?<=[.!?])\s+', clean)[:150]:
        if len(s)>5: db.session.add(Sentence(en=s, pt="...", story_id=sid))
    db.session.commit()
    return redirect(url_for('ler', id=sid))

# --- ROTAS PADRÃO ---
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

@app.route('/search', methods=['GET', 'POST'])
def search():
    res = []
    if request.method == 'POST': res = utils.search_pubmed(request.form.get('query'))
    return render_template('layout.html', mode='search', results=res, app_name=APP_NAME)

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

@app.route('/pronunciation')
def pronunciation():
    c = Card.query.order_by(db.func.random()).limit(10).all()
    w = [{"w": x.front, "ipa": x.ipa} for x in c] if c else [{"w":"Science", "ipa":"/ˈsaɪ.əns/"}]
    return render_template('layout.html', mode='pronunciation', words_json=json.dumps(w), app_name=APP_NAME)

@app.route('/writer')
def writer(): return render_template('layout.html', mode='writer', connectors=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

@app.route('/checker', methods=['GET', 'POST'])
def checker():
    o, c = "", ""
    if request.method == 'POST': o=request.form.get('text_input'); c=utils.improve_english_text(o)
    return render_template('layout.html', mode='checker', original=o, corrected=c, app_name=APP_NAME)

@app.route('/miner', methods=['GET', 'POST'])
def miner():
    k = []
    if request.method == 'POST':
        t = request.form.get('text_input', '')
        f = request.files.get('file_input')
        if f: 
            try: t += " ".join([p.extract_text() for p in PdfReader(f).pages])
            except: pass
        k = utils.get_top_keywords(t)
    return render_template('layout.html', mode='miner', keywords=k, app_name=APP_NAME)

@app.route('/summarizer', methods=['GET', 'POST'])
def summarizer():
    if request.method == 'POST':
        t = request.form.get('text_input', '')
        f = request.files.get('file_input')
        if f: 
            try: t = " ".join([p.extract_text() for p in PdfReader(f).pages])
            except: pass
        sid = str(uuid.uuid4())[:8]
        db.session.add(Story(id=sid, title="Resumo: "+t[:20]))
        clean = re.sub(r'<[^>]*>', '', utils.format_abstract_smart(t))
        for s in re.split(r'(?<=[.!?])\s+', clean)[:50]:
            if len(s)>5: db.session.add(Sentence(en=s, pt="...", story_id=sid))
        db.session.commit()
        return redirect(url_for('ler', id=sid))
    return render_template('layout.html', mode='summarizer', app_name=APP_NAME)

@app.route('/delete_story/<id>')
def delete_story(id): db.session.delete(Story.query.get(id)); db.session.commit(); return redirect(url_for('index'))
@app.route('/delete_ref/<int:id>')
def delete_ref(id): db.session.delete(Reference.query.get(id)); db.session.commit(); return redirect(url_for('library'))

if __name__ == '__main__': app.run(host='0.0.0.0', port=5000)
