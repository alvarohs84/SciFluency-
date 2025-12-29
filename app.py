import os
import uuid
import json
from datetime import datetime, timedelta
from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from werkzeug.utils import secure_filename
from sqlalchemy import text, inspect
from pypdf import PdfReader
from deep_translator import GoogleTranslator

# Importa modelos e utils
from models import db, Deck, Card, Story, Sentence, StudyLog
import utils

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = os.environ.get('DATABASE_URL', 'sqlite:///scifluency.db')
if app.config['SQLALCHEMY_DATABASE_URI'].startswith("postgres://"):
    app.config['SQLALCHEMY_DATABASE_URI'] = app.config['SQLALCHEMY_DATABASE_URI'].replace("postgres://", "postgresql://", 1)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

db.init_app(app)
APP_NAME = "SciFluency"

# --- AUTO-REPARO (Para garantir que o banco não quebre) ---
def fix_db():
    with app.app_context():
        try:
            db.create_all()
            # Garante tabelas básicas
            if not Deck.query.get("my_vocab"):
                db.session.add(Deck(id="my_vocab", name="Meu Vocabulário", icon="🗂️"))
                db.session.commit()
        except: pass
fix_db()

@app.route('/')
def index():
    # Tela inicial limpa: Mostra textos salvos e estatísticas de palavras
    stories = Story.query.order_by(Story.id.desc()).all()
    vocab_count = Card.query.count()
    return render_template('layout.html', mode='home', stories=stories, vocab_count=vocab_count, app_name=APP_NAME)

@app.route('/novo', methods=['GET', 'POST'])
def novo():
    if request.method == 'POST':
        # Processador Universal (PDF ou Texto)
        full_text = request.form.get('text_input', '')
        f = request.files.get('file_input')
        
        if f and f.filename:
            if f.filename.endswith('.pdf'):
                try: full_text = " ".join([p.extract_text() for p in PdfReader(f).pages])
                except: pass
        
        if full_text.strip():
            sid = str(uuid.uuid4())[:8]
            title = full_text[:50] + "..."
            db.session.add(Story(id=sid, title=title))
            
            # Formatação inteligente é feita no layout agora
            # Salvamos apenas a frase bruta no banco para simplicidade
            import re
            sentences = re.split(r'(?<=[.!?])\s+', full_text)
            tr = GoogleTranslator(source='en', target='pt')
            
            for s in sentences[:100]: # Limite para não travar
                if len(s) > 3:
                    try: pt = tr.translate(s)
                    except: pt = "..."
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
    # Modo Focado em Fala
    cards = Card.query.order_by(db.func.random()).limit(15).all()
    words = [{"w": c.front, "ipa": c.ipa} for c in cards]
    if not words: words = [{"w": "Welcome", "ipa": "/ˈwɛlkəm/"}]
    return render_template('layout.html', mode='pronunciation', words_json=json.dumps(words), app_name=APP_NAME)

@app.route('/vocab')
def vocab():
    cards = Card.query.order_by(Card.id.desc()).all()
    return render_template('layout.html', mode='vocab', cards=cards, app_name=APP_NAME)

@app.route('/tts', methods=['POST'])
def tts():
    audio = utils.get_audio_sync(request.form.get('text'), request.form.get('accent'))
    if audio: return send_file(audio, mimetype='audio/mpeg')
    return "Error", 500

@app.route('/translate')
def translate():
    w = request.args.get('w')
    t = GoogleTranslator(source='en', target='pt').translate(w)
    return jsonify({'t': t})

@app.route('/add_card', methods=['POST'])
def add_card():
    front = request.form.get('front')
    back = request.form.get('back')
    ipa_txt = utils.get_phonetic(front)
    db.session.add(Card(front=front, back=back, ipa=ipa_txt, deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d')))
    db.session.commit()
    return jsonify({'status': 'ok'})

@app.route('/delete_story/<id>')
def delete_story(id):
    s = Story.query.get(id)
    db.session.delete(s); db.session.commit()
    return redirect(url_for('index'))

@app.route('/delete_card/<int:id>')
def delete_card(id):
    c = Card.query.get(id)
    db.session.delete(c); db.session.commit()
    return redirect(url_for('vocab'))

if __name__ == '__main__': app.run(host='0.0.0.0', port=5000)
