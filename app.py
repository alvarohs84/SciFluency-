import os
import re
from datetime import datetime, timedelta
from io import BytesIO

from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from werkzeug.utils import secure_filename
from sqlalchemy import text, inspect  # Adicionado inspect
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
APP_NAME = "SciFluency"

# --- AUTO-REPARO DO BANCO DE DADOS (ROBUSTO) ---
def fix_database_schema():
    """Verifica e corrige colunas faltantes no banco de dados com Inspect"""
    with app.app_context():
        try:
            # Garante que as tabelas existem
            db.create_all()
            
            inspector = inspect(db.engine)
            
            # --- CORREÇÃO TABELA CARD ---
            if inspector.has_table('card'):
                columns = [c['name'] for c in inspector.get_columns('card')]
                with db.engine.connect() as conn:
                    # Adiciona ease_factor se faltar
                    if 'ease_factor' not in columns:
                        print("🛠️ Migrando: Adicionando ease_factor...")
                        conn.execute(text("ALTER TABLE card ADD COLUMN ease_factor FLOAT DEFAULT 2.5"))
                        conn.commit()
                    
                    # Adiciona context se faltar (para garantir)
                    if 'context' not in columns:
                        print("🛠️ Migrando: Adicionando context...")
                        conn.execute(text("ALTER TABLE card ADD COLUMN context TEXT"))
                        conn.commit()

            print("✅ Banco de dados verificado e atualizado.")
            
        except Exception as e:
            print(f"❌ Erro crítico na migração: {e}")

# Executa o reparo imediatamente ao iniciar o app
fix_database_schema()

# --- ROTAS PRINCIPAIS ---

@app.route('/')
def index():
    # Inicialização segura de dados padrão
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
    
    # Estatísticas
    try:
        stats = {
            "refs": Reference.query.count(),
            "to_read": Reference.query.filter_by(status='to_read').count(),
            "vocab": Card.query.count(),
            "mastered": Card.query.filter(Card.interval > 21).count(),
            "heatmap": []
        }
        # Heatmap seguro
        for i in range(6, -1, -1):
            d = (datetime.now() - timedelta(days=i)).strftime('%Y-%m-%d')
            l = StudyLog.query.filter_by(date=d).first()
            stats['heatmap'].append({"date": d, "color": "var(--learn)" if l else "#ecf0f1"})
    except:
        stats = {"refs":0, "to_read":0, "vocab":0, "mastered":0, "heatmap":[]}

    drafts = Draft.query.filter_by(project_id=proj.id).all()
    return render_template('layout.html', mode='dashboard', project=proj, drafts=drafts, stats=stats, app_name=APP_NAME)

@app.route('/library', methods=['GET', 'POST'])
def library():
    if request.method == 'POST':
        f = request.files.get('pdf_file')
        if f and f.filename:
            filename = secure_filename(f.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            f.save(filepath)
            try:
                reader = PdfReader(filepath)
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
    ref = Reference.query.get(id)
    if not ref: return redirect(url_for('library'))
    try:
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], ref.pdf_filename)
        reader = PdfReader(filepath)
        full_text = " ".join([page.extract_text() or "" for page in reader.pages])
        formatted_html = utils.format_abstract_smart(full_text)
        if ref.status == 'to_read':
            ref.status = 'reading'
            db.session.commit()
        return render_template('layout.html', mode='reader_tool', title=ref.title, content=formatted_html, app_name=APP_NAME)
    except Exception as e: return f"Erro ao ler PDF: {e}"

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
    today = datetime.now().strftime('%Y-%m-%d')
    log = StudyLog.query.filter_by(date=today).first()
    if not log: log = StudyLog(date=today, count=0); db.session.add(log)
    log.count += 1
    db.session.commit()
    return jsonify({"status": "ok"})

@app.route('/jogar/<id>')
def jogar(id):
    target_deck = "my_vocab"
    # Lógica segura para SRS
    due = Card.query.filter(Card.deck_id == target_deck, Card.next_review <= datetime.now().strftime('%Y-%m-%d')).all()
    card_data = None
    if due:
        c = due[0]
        card_data = {"id": c.id, "front": c.front, "back": c.back, "ipa": c.ipa}
    return render_template('layout.html', mode='study_play', deck=Deck.query.get(target_deck), card=card_data, app_name=APP_NAME)

@app.route('/rate_card', methods=['POST'])
def rate_card():
    c = Card.query.filter_by(deck_id="my_vocab", front=request.form.get('front')).first()
    if c:
        rating = request.form.get('rating')
        if rating == 'hard': c.interval = 1
        elif rating == 'good': c.interval = max(1, int(c.interval * 1.5))
        elif rating == 'easy': c.interval = max(2, int(c.interval * 2.5))
        c.next_review = (datetime.now() + timedelta(days=c.interval)).strftime('%Y-%m-%d')
        db.session.commit()
    return redirect(url_for('jogar', id="my_vocab"))

# --- FERRAMENTAS EXTRAS ---
@app.route('/tutor', methods=['GET', 'POST'])
def tutor():
    res = None
    if request.method == 'POST':
        pt = request.form.get('pt_text')
        en = GoogleTranslator(source='pt', target='en').translate(pt)
        res = {"en": en, "pt": pt}
    return render_template('layout.html', mode='tutor', res=res, app_name=APP_NAME)

@app.route('/pronunciation')
def pronunciation():
    # Palavras de exemplo se não houver cards
    words = [{"w": c.front, "ipa": c.ipa} for c in Card.query.limit(10).all()]
    if not words: words = [{"w": "Research", "ipa": "/rɪˈsɜːrtʃ/"}]
    return render_template('layout.html', mode='pronunciation', words_json=utils.json.dumps(words), app_name=APP_NAME)

@app.route('/phrases')
def phrases(): return render_template('layout.html', mode='phrases', phrases=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

@app.route('/novo')
def novo(): return render_template('layout.html', mode='tutor', res=None, app_name=APP_NAME) # Reusa layout tutor para input manual

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
