import os
import re
from datetime import datetime, timedelta
from io import BytesIO

from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from werkzeug.utils import secure_filename
from sqlalchemy import text
from pypdf import PdfReader
from deep_translator import GoogleTranslator
from gtts import gTTS
import edge_tts
import asyncio

# IMPORTAÇÕES LOCAIS (Models e Utils mantidos)
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
APP_NAME = "SciFluency"

# --- ROTAS PRINCIPAIS ---

@app.route('/')
def index():
    # Garante que o usuário tenha um projeto e um deck inicial
    if not Project.query.first():
        try:
            db.session.add(Project(id="thesis", title="Meu Projeto Acadêmico", target_journal="A Definir"))
            sections = ["1. Introduction", "2. Literature Review", "3. Methods", "4. Results", "5. Discussion", "6. Conclusion"]
            for sec in sections:
                db.session.add(Draft(section_name=sec, content="", project_id="thesis"))
            if not Deck.query.get("my_vocab"):
                db.session.add(Deck(id="my_vocab", name="Vocabulário Acadêmico", icon="🎓"))
            db.session.commit()
        except: db.session.rollback()
    
    proj = Project.query.first()
    
    # Estatísticas Unificadas
    stats = {
        "refs": Reference.query.count(),
        "to_read": Reference.query.filter_by(status='to_read').count(),
        "vocab": Card.query.count(),
        "mastered": Card.query.filter(Card.interval > 21).count(),
        "heatmap": [] # Heatmap simplificado
    }
    
    # Gera heatmap dos últimos 7 dias
    for i in range(6, -1, -1):
        d = (datetime.now() - timedelta(days=i)).strftime('%Y-%m-%d')
        l = StudyLog.query.filter_by(date=d).first()
        stats['heatmap'].append({"date": d, "color": "var(--learn)" if l else "#ecf0f1"})

    drafts = Draft.query.filter_by(project_id=proj.id).all()
    return render_template('layout.html', mode='dashboard', project=proj, drafts=drafts, stats=stats, app_name=APP_NAME)

# --- BIBLIOTECA E LEITURA DE PDFS ---

@app.route('/library', methods=['GET', 'POST'])
def library():
    if request.method == 'POST':
        f = request.files.get('pdf_file')
        if f and f.filename:
            filename = secure_filename(f.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            f.save(filepath)
            
            # Extrai texto para resumo inicial
            try:
                reader = PdfReader(filepath)
                # Tenta pegar apenas a primeira página para o resumo
                text_content = reader.pages[0].extract_text()
                abstract_preview = text_content[:400] + "..."
            except: abstract_preview = "Texto não detectado."

            new_ref = Reference(
                title=request.form.get('title') or filename,
                authors=request.form.get('authors') or "Unknown",
                year=request.form.get('year') or "2024",
                status='to_read', pdf_filename=filename, abstract=abstract_preview, project_id="thesis"
            )
            db.session.add(new_ref); db.session.commit()
            return redirect(url_for('library'))

    refs = Reference.query.order_by(Reference.id.desc()).all()
    return render_template('layout.html', mode='library', references=refs, app_name=APP_NAME)

@app.route('/read_ref/<int:id>')
def read_ref(id):
    """Lê um PDF da biblioteca com as ferramentas de áudio e dicionário"""
    ref = Reference.query.get(id)
    if not ref: return redirect(url_for('library'))
    
    try:
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], ref.pdf_filename)
        reader = PdfReader(filepath)
        full_text = " ".join([page.extract_text() or "" for page in reader.pages])
        
        # Usa a função inteligente para formatar o texto para leitura (limpa quebras, bold em seções)
        formatted_html = utils.format_abstract_smart(full_text)
        
        # Atualiza status para 'lendo'
        if ref.status == 'to_read':
            ref.status = 'reading'
            db.session.commit()
            
        return render_template('layout.html', mode='reader_tool', title=ref.title, content=formatted_html, app_name=APP_NAME)
        
    except Exception as e:
        return f"Erro ao ler PDF: {e}"

@app.route('/update_status/<int:id>/<new_status>')
def update_status(id, new_status):
    ref = Reference.query.get(id)
    if ref: ref.status = new_status; db.session.commit()
    return redirect(url_for('library'))

@app.route('/delete_ref/<int:id>')
def delete_ref(id):
    ref = Reference.query.get(id)
    if ref:
        try: os.remove(os.path.join(app.config['UPLOAD_FOLDER'], ref.pdf_filename))
        except: pass
        db.session.delete(ref); db.session.commit()
    return redirect(url_for('library'))

# --- OUTROS (Mantidos Integralmente) ---

@app.route('/writer/<int:draft_id>', methods=['GET', 'POST'])
def writer(draft_id):
    draft = Draft.query.get(draft_id)
    if request.method == 'POST':
        draft.content = request.form.get('content')
        draft.last_updated = datetime.utcnow()
        db.session.commit()
        return jsonify({"status": "saved", "time": datetime.now().strftime("%H:%M")})
    return render_template('layout.html', mode='writer', draft=draft, connectors=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

@app.route('/vocab_list')
def vocab_list(): return render_template('layout.html', mode='vocab_list', cards=Card.query.all(), app_name=APP_NAME)

@app.route('/search', methods=['GET', 'POST'])
def search(): return render_template('layout.html', mode='search', results=[], links=utils.RESEARCH_LINKS, app_name=APP_NAME)

@app.route('/tts', methods=['POST'])
def tts():
    audio = utils.get_audio_sync(request.form.get('text'), request.form.get('accent'))
    return send_file(audio, mimetype='audio/mpeg') if audio else ("Error", 500)

@app.route('/traduzir_palavra')
def traduzir(): return jsonify({"t": GoogleTranslator(source='en', target='pt').translate(request.args.get('w'))})

@app.route('/adicionar_vocab', methods=['POST'])
def add_vocab():
    db.session.add(Card(front=request.form.get('term'), back="...", deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d')))
    
    # Loga atividade (Gamification)
    today = datetime.now().strftime('%Y-%m-%d')
    log = StudyLog.query.filter_by(date=today).first()
    if not log: log = StudyLog(date=today, count=0); db.session.add(log)
    log.count += 1
    
    db.session.commit()
    return jsonify({"status": "ok"})

@app.route('/jogar/<id>')
def jogar(id):
    target_deck = "my_vocab"
    due = Card.query.filter(Card.deck_id == target_deck, Card.next_review <= datetime.now().strftime('%Y-%m-%d')).all()
    card_data = None
    if due:
        c = due[0] # Pega o primeiro da fila
        card_data = {"id": c.id, "front": c.front, "back": c.back, "ipa": c.ipa}
    return render_template('layout.html', mode='study_play', deck=Deck.query.get(target_deck), card=card_data, app_name=APP_NAME)

@app.route('/rate_card', methods=['POST'])
def rate_card():
    # Lógica SRS Simplificada
    c = Card.query.filter_by(deck_id="my_vocab", front=request.form.get('front')).first()
    if c:
        rating = request.form.get('rating')
        if rating == 'hard': c.interval = 1
        elif rating == 'good': c.interval = max(1, int(c.interval * 1.5))
        elif rating == 'easy': c.interval = max(2, int(c.interval * 2.5))
        c.next_review = (datetime.now() + timedelta(days=c.interval)).strftime('%Y-%m-%d')
        db.session.commit()
    return redirect(url_for('jogar', id="my_vocab"))

# --- AUTO-REPARO ---
def fix_database_schema():
    with app.app_context():
        try:
            db.create_all()
            with db.engine.connect() as conn:
                try: conn.execute(text("SELECT ease_factor FROM card LIMIT 1"))
                except: 
                    try: conn.execute(text("ALTER TABLE card ADD COLUMN ease_factor FLOAT DEFAULT 2.5")); conn.commit()
                    except: pass
        except: pass

fix_database_schema()

if __name__ == '__main__': app.run(host='0.0.0.0', port=5000)
