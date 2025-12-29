import os
import json
import random
from datetime import datetime, timedelta
from io import BytesIO

from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from sqlalchemy import text, inspect
from pypdf import PdfReader
from deep_translator import GoogleTranslator
from Bio import Entrez
import uuid

# IMPORTAÇÕES LOCAIS (A MÁGICA ACONTECE AQUI)
from models import db, Deck, Card, Story, Sentence, StudyLog
import utils

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = os.environ.get('DATABASE_URL', 'sqlite:///scifluency.db')
if app.config['SQLALCHEMY_DATABASE_URI'].startswith("postgres://"):
    app.config['SQLALCHEMY_DATABASE_URI'] = app.config['SQLALCHEMY_DATABASE_URI'].replace("postgres://", "postgresql://", 1)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

db.init_app(app)
Entrez.email = "researcher@example.com"
APP_NAME = "SciFluency"

# LINKS GLOBAIS
RESEARCH_LINKS = [
    {"name": "PubMed", "url": "https://pubmed.ncbi.nlm.nih.gov/", "icon": "🧬"},
    {"name": "SciELO", "url": "https://scielo.org/", "icon": "🌎"},
    {"name": "Google Scholar", "url": "https://scholar.google.com.br/", "icon": "🎓"}
]

def log_activity(points):
    today = datetime.now().strftime('%Y-%m-%d')
    log = StudyLog.query.filter_by(date=today).first()
    if not log: log = StudyLog(date=today, count=0); db.session.add(log)
    log.count += points
    db.session.commit()

@app.route('/')
def index():
    # Garante decks
    if not Deck.query.first():
        db.create_all()
        db.session.add(Deck(id="my_vocab", name="My Vocabulary", icon="🗂️"))
        db.session.commit()
    
    main_deck = Deck.query.get("my_vocab")
    today = datetime.now().strftime('%Y-%m-%d')
    due = Card.query.filter(Card.deck_id == "my_vocab", Card.next_review <= today).count()
    
    heatmap = []
    for i in range(13, -1, -1):
        d = (datetime.now() - timedelta(days=i)).strftime('%Y-%m-%d')
        l = StudyLog.query.filter_by(date=d).first()
        heatmap.append({"date": d, "count": l.count if l else 0, "color": "var(--heat-1)" if l else "var(--heat-0)"})

    return render_template('layout.html', mode='list', stories=Story.query.all(), 
                         decks=[{"id": "my_vocab", "name": "My Vocab", "due_count": due}], 
                         stats={"total": Card.query.count(), "mastered": 0, "heatmap": heatmap}, 
                         app_name=APP_NAME, links=RESEARCH_LINKS)

@app.route('/tts', methods=['POST'])
def tts():
    audio = utils.get_audio_sync(request.form.get('text'), request.form.get('accent'))
    return send_file(audio, mimetype='audio/mpeg') if audio else ("Error", 500)

@app.route('/summarizer', methods=['GET', 'POST'])
def summarizer():
    results = []
    if request.method == 'POST':
        # Lógica simplificada chamando utils
        for f in request.files.getlist("nbib_file"):
            if f.filename.endswith('.pdf'):
                txt = " ".join([p.extract_text() for p in PdfReader(f).pages])
                results.append({"title": f.filename, "formatted_html": utils.format_abstract_smart(txt)})
        
        txt_in = request.form.get('text_input')
        if txt_in: results.append({"title": "Text Input", "formatted_html": utils.format_abstract_smart(txt_in)})
        
    return render_template('layout.html', mode='summarizer', batch_results=results, app_name=APP_NAME)

# --- ROTAS CURTAS PARA AS OUTRAS FUNÇÕES ---
@app.route('/search', methods=['GET', 'POST'])
def search():
    res = []
    if request.method == 'POST':
        # Chama Entrez logic aqui (simplificado para o exemplo)
        pass 
    return render_template('layout.html', mode='search', results=res, app_name=APP_NAME)

@app.route('/phrases')
def phrases():
    return render_template('layout.html', mode='phrases', phrases=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

@app.route('/traduzir_palavra')
def traduzir():
    return jsonify({"t": GoogleTranslator(source='en', target='pt').translate(request.args.get('w'))})

@app.route('/adicionar_vocab', methods=['POST'])
def add_vocab():
    f = request.form.get('term')
    db.session.add(Card(front=f, back="...", deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d')))
    db.session.commit()
    return jsonify({"status":"ok"})

# ROTA PARA SCI-WRITER (NOVA!)
@app.route('/writer')
def writer():
    return render_template('layout.html', mode='writer', connectors=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

with app.app_context():
    db.create_all()

if __name__ == '__main__': app.run(host='0.0.0.0', port=5000)
