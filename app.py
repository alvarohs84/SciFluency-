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

# IMPORTAÇÕES LOCAIS
from models import db, Deck, Card, Story, Sentence, StudyLog, Reference, Project, Note, Draft
import utils

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = os.environ.get('DATABASE_URL', 'sqlite:///scifluency.db')
if app.config['SQLALCHEMY_DATABASE_URI'].startswith("postgres://"):
    app.config['SQLALCHEMY_DATABASE_URI'] = app.config['SQLALCHEMY_DATABASE_URI'].replace("postgres://", "postgresql://", 1)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# CONFIGURAÇÃO DE UPLOAD
UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

db.init_app(app)
APP_NAME = "SciFluency Research OS"

# --- ROTAS PRINCIPAIS ---

@app.route('/')
def index():
    # Cria projeto padrão se não existir
    if not Project.query.first():
        try:
            db.session.add(Project(id="thesis", title="Meu Projeto de Pesquisa", target_journal="Nature/Science"))
            # Verifica se o deck existe antes de criar
            if not Deck.query.get("my_vocab"):
                db.session.add(Deck(id="my_vocab", name="Vocabulário Acadêmico", icon="🎓"))
            db.session.commit()
        except:
            db.session.rollback()
    
    proj = Project.query.first()
    # Contagem de Referências
    refs_count = Reference.query.count()
    to_read = Reference.query.filter_by(status='to_read').count()
    reading = Reference.query.filter_by(status='reading').count()
    done = Reference.query.filter_by(status='done').count()
    
    stats = {
        "refs": refs_count, "to_read": to_read, "reading": reading, "done": done,
        "vocab": Card.query.count()
    }

    return render_template('layout.html', mode='dashboard', project=proj, stats=stats, app_name=APP_NAME)

@app.route('/library', methods=['GET', 'POST'])
def library():
    if request.method == 'POST':
        # UPLOAD DE REFERÊNCIA (PDF)
        f = request.files.get('pdf_file')
        if f and f.filename:
            filename = secure_filename(f.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            f.save(filepath)
            
            # Tenta extrair texto para o resumo automático
            try:
                reader = PdfReader(filepath)
                text_content = " ".join([page.extract_text() for page in reader.pages[:2]]) # Pega as 2 primeiras pags
                abstract_preview = text_content[:500] + "..."
            except:
                abstract_preview = "Texto não extraído."

            new_ref = Reference(
                title=request.form.get('title') or filename,
                authors=request.form.get('authors') or "Unknown",
                year=request.form.get('year') or "2024",
                status='to_read',
                pdf_filename=filename,
                abstract=abstract_preview,
                project_id="thesis"
            )
            db.session.add(new_ref)
            db.session.commit()
            return redirect(url_for('library'))

    refs = Reference.query.order_by(Reference.id.desc()).all()
    return render_template('layout.html', mode='library', references=refs, app_name=APP_NAME)

@app.route('/update_status/<int:id>/<new_status>')
def update_status(id, new_status):
    ref = Reference.query.get(id)
    if ref:
        ref.status = new_status
        db.session.commit()
    return redirect(url_for('library'))

@app.route('/delete_ref/<int:id>')
def delete_ref(id):
    ref = Reference.query.get(id)
    if ref:
        # Tenta remover o arquivo físico também
        try: os.remove(os.path.join(app.config['UPLOAD_FOLDER'], ref.pdf_filename))
        except: pass
        db.session.delete(ref)
        db.session.commit()
    return redirect(url_for('library'))

# --- RECURSOS DE ESTUDO (MANTIDOS) ---
@app.route('/vocab_list')
def vocab_list():
    return render_template('layout.html', mode='vocab_list', cards=Card.query.all(), app_name=APP_NAME)

@app.route('/writer')
def writer():
    return render_template('layout.html', mode='writer', connectors=utils.ACADEMIC_PHRASEBANK, app_name=APP_NAME)

@app.route('/search', methods=['GET', 'POST'])
def search():
    results = []
    if request.method == 'POST':
        # Simulação de busca no PubMed (pode ser reativada completa com BioPython)
        pass 
    return render_template('layout.html', mode='search', results=results, links=utils.RESEARCH_LINKS, app_name=APP_NAME)

@app.route('/tts', methods=['POST'])
def tts():
    # Rota para voz neural
    text = request.form.get('text')
    accent = request.form.get('accent')
    audio = utils.get_audio_sync(text, accent)
    return send_file(audio, mimetype='audio/mpeg') if audio else ("Error", 500)

@app.route('/traduzir_palavra')
def traduzir():
    return jsonify({"t": GoogleTranslator(source='en', target='pt').translate(request.args.get('w'))})

@app.route('/adicionar_vocab', methods=['POST'])
def add_vocab():
    term = request.form.get('term')
    db.session.add(Card(front=term, back="...", deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d')))
    db.session.commit()
    return jsonify({"status": "ok"})

# --- AUTO-REPARO DO BANCO DE DADOS ---
def fix_database_schema():
    """Verifica e corrige colunas faltantes no banco de dados"""
    with app.app_context():
        try:
            # 1. Cria tabelas novas (Project, Reference, etc)
            db.create_all()
            
            # 2. Força a adição da coluna 'ease_factor' se ela não existir
            # Isso resolve o erro 'UndefinedColumn'
            with db.engine.connect() as conn:
                # Tenta selecionar a coluna para ver se existe
                try:
                    conn.execute(text("SELECT ease_factor FROM card LIMIT 1"))
                except:
                    # Se der erro, é porque não existe. Adiciona.
                    # A sintaxe varia entre SQLite e Postgres, este comando tenta ser compatível ou foca no Postgres
                    try:
                        conn.execute(text("ALTER TABLE card ADD COLUMN ease_factor FLOAT DEFAULT 2.5"))
                        conn.commit()
                        print("✅ Coluna 'ease_factor' adicionada com sucesso.")
                    except Exception as e:
                        print(f"⚠️ Aviso na migração: {e}")
        except Exception as e:
            print(f"Erro geral no banco: {e}")

# Executa o reparo ao iniciar
fix_database_schema()

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
