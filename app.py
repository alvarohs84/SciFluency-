import os
import json
import csv
import uuid
import re
import random
import time
from io import BytesIO, StringIO, TextIOWrapper
from datetime import datetime, timedelta

from flask import Flask, render_template, request, redirect, url_for, send_file, jsonify
from sqlalchemy import text, inspect
from gtts import gTTS
from Bio import Entrez
from deep_translator import GoogleTranslator
from pypdf import PdfReader

# --- IMPORTAÇÕES LOCAIS ---
from models import db, Deck, Card, Story, Sentence, StudyLog
import utils

app = Flask(__name__)

# Configuração do Banco de Dados
app.config['SQLALCHEMY_DATABASE_URI'] = os.environ.get('DATABASE_URL', 'sqlite:///scifluency.db')
if app.config['SQLALCHEMY_DATABASE_URI'].startswith("postgres://"):
    app.config['SQLALCHEMY_DATABASE_URI'] = app.config['SQLALCHEMY_DATABASE_URI'].replace("postgres://", "postgresql://", 1)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

Entrez.email = "researcher@example.com"
APP_NAME = "SciFluency"

db.init_app(app)

# --- LINKS DO RESEARCH HUB ---
RESEARCH_LINKS = [
    {"name": "PubMed", "url": "https://pubmed.ncbi.nlm.nih.gov/", "icon": "🧬", "desc": "Biomedical Literature"},
    {"name": "SciELO", "url": "https://scielo.org/", "icon": "🌎", "desc": "Open Access Journals"},
    {"name": "Google Scholar", "url": "https://scholar.google.com.br/", "icon": "🎓", "desc": "Academic Search"},
    {"name": "Cochrane", "url": "https://www.cochranelibrary.com/", "icon": "🏥", "desc": "Evidence-Based Medicine"},
    {"name": "Scopus", "url": "https://www.scopus.com/", "icon": "🔭", "desc": "Citation Database"},
    {"name": "Web of Science", "url": "https://www.webofscience.com/", "icon": "🕸️", "desc": "Scientific Citation Index"}
]

# --- BANCO DE FRASES ACADÊMICAS (O RECHEIO) ---
ACADEMIC_PHRASEBANK = {
    "1. Introduction & Context": [
        {"en": "Recent developments in this field have heightened the need for...", "pt": "Desenvolvimentos recentes neste campo aumentaram a necessidade de..."},
        {"en": "Currently, there is a paucity of data regarding...", "pt": "Atualmente, há escassez de dados sobre..."},
        {"en": "This study aims to investigate the relationship between...", "pt": "Este estudo visa investigar a relação entre..."},
        {"en": "Previous research has established that...", "pt": "Pesquisas anteriores estabeleceram que..."},
        {"en": "However, these results remain controversial.", "pt": "No entanto, esses resultados permanecem controversos."},
        {"en": "The primary objective of this paper is to evaluate...", "pt": "O objetivo principal deste artigo é avaliar..."},
        {"en": "Understanding the mechanisms underlying this phenomenon is crucial.", "pt": "Compreender os mecanismos subjacentes a este fenômeno é crucial."},
        {"en": "Little is known about the effects of...", "pt": "Pouco se sabe sobre os efeitos de..."}
    ],
    "2. Methods & Materials": [
        {"en": "Data were collected using a semi-structured interview guide.", "pt": "Os dados foram coletados usando um roteiro de entrevista semiestruturado."},
        {"en": "The participants were divided into two groups.", "pt": "Os participantes foram divididos em dois grupos."},
        {"en": "Statistical analysis was performed using SPSS software.", "pt": "A análise estatística foi realizada usando o software SPSS."},
        {"en": "We excluded patients with a history of...", "pt": "Excluímos pacientes com histórico de..."},
        {"en": "Samples were obtained from...", "pt": "As amostras foram obtidas de..."},
        {"en": "The experiment was conducted in triplicate.", "pt": "O experimento foi conduzido em triplicata."},
        {"en": "To adjust for potential confounders, we used...", "pt": "Para ajustar potenciais confundidores, usamos..."},
        {"en": "All subjects provided informed consent.", "pt": "Todos os sujeitos forneceram consentimento informado."}
    ],
    "3. Results & Findings": [
        {"en": "There was a significant correlation between...", "pt": "Houve uma correlação significativa entre..."},
        {"en": "Table 1 presents the demographic characteristics of the sample.", "pt": "A Tabela 1 apresenta as características demográficas da amostra."},
        {"en": "Contrary to expectations, we did not find...", "pt": "Contrário às expectativas, não encontramos..."},
        {"en": "The results indicate that...", "pt": "Os resultados indicam que..."},
        {"en": "Interestingly, no significant difference was observed.", "pt": "Curiosamente, nenhuma diferença significativa foi observada."},
        {"en": "Figure 2 illustrates the breakdown of...", "pt": "A Figura 2 ilustra a distribuição de..."},
        {"en": "Our findings are consistent with those of previous studies.", "pt": "Nossos achados são consistentes com os de estudos anteriores."},
        {"en": "A clear trend was identified regarding...", "pt": "Uma tendência clara foi identificada em relação a..."}
    ],
    "4. Discussion & Argumentation": [
        {"en": "These findings suggest that...", "pt": "Esses achados sugerem que..."},
        {"en": "One possible explanation for this result is...", "pt": "Uma possível explicação para este resultado é..."},
        {"en": "This study provides new insights into...", "pt": "Este estudo fornece novos insights sobre..."},
        {"en": "However, some limitations should be noted.", "pt": "No entanto, algumas limitações devem ser notadas."},
        {"en": "Our results challenge the widely held view that...", "pt": "Nossos resultados desafiam a visão amplamente aceita de que..."},
        {"en": "It is plausible to assume that...", "pt": "É plausível assumir que..."},
        {"en": "Further research is needed to confirm these findings.", "pt": "Mais pesquisas são necessárias para confirmar esses achados."},
        {"en": "The implications of this study are twofold.", "pt": "As implicações deste estudo são duplas."}
    ],
    "5. Conclusion": [
        {"en": "In conclusion, this study demonstrates that...", "pt": "Em conclusão, este estudo demonstra que..."},
        {"en": "The evidence from this study suggests...", "pt": "As evidências deste estudo sugerem..."},
        {"en": "Ideally, future studies should address...", "pt": "Idealmente, estudos futuros devem abordar..."},
        {"en": "Overall, these results highlight the importance of...", "pt": "No geral, esses resultados destacam a importância de..."},
        {"en": "We recommend that...", "pt": "Recomendamos que..."}
    ]
}

# --- FUNÇÕES DE APOIO ---
def log_activity(points):
    today = datetime.now().strftime('%Y-%m-%d')
    log = StudyLog.query.filter_by(date=today).first()
    if not log: log = StudyLog(date=today, count=0); db.session.add(log)
    log.count += points
    db.session.commit()

def seed_database():
    if Deck.query.filter_by(id="my_vocab").first() is None:
        db.session.add(Deck(id="my_vocab", name="My Vocabulary", pt_name="Meu Banco de Palavras", icon="🗂️"))
        db.session.add(Card(front="Welcome", back="Bem-vindo", ipa="/ˈwel.kəm/", context="Welcome to SciFluency.", deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d')))
        db.session.commit()

def check_and_migrate_db():
    try:
        inspector = inspect(db.engine)
        if inspector.has_table('card'):
            columns = [c['name'] for c in inspector.get_columns('card')]
            if 'context' not in columns:
                with db.engine.connect() as conn:
                    conn.execute(text("ALTER TABLE card ADD COLUMN context TEXT"))
                    conn.commit()
    except Exception as e: print(f"Migration warning: {e}")

# --- ROTAS ---
@app.route('/')
def index():
    main_deck = Deck.query.get("my_vocab")
    if not main_deck: seed_database(); main_deck = Deck.query.get("my_vocab")
    today = datetime.now().strftime('%Y-%m-%d')
    due = Card.query.filter(Card.deck_id == "my_vocab", Card.next_review <= today).count()
    deck_data = {"id": main_deck.id, "name": main_deck.name, "pt": main_deck.pt_name, "icon": main_deck.icon, "due_count": due}
    
    total_cards = Card.query.count()
    mastered_cards = Card.query.filter(Card.interval > 21).count()
    heatmap = []
    for i in range(13, -1, -1):
        d_date = (datetime.now() - timedelta(days=i)).strftime('%Y-%m-%d')
        log = StudyLog.query.filter_by(date=d_date).first()
        count = log.count if log else 0
        color = "var(--heat-0)"
        if count > 0: color = "var(--heat-1)"
        if count > 5: color = "var(--heat-2)"
        if count > 15: color = "var(--heat-3)"
        heatmap.append({"date": d_date, "count": count, "color": color, "day": d_date[-2:]})
    stats = {"total": total_cards, "mastered": mastered_cards, "heatmap": heatmap}
    return render_template('layout.html', mode='list', stories=Story.query.all(), decks=[deck_data], stats=stats, app_name=APP_NAME, links=RESEARCH_LINKS)

@app.route('/reset_decks')
def reset_decks():
    db.drop_all(); db.create_all(); seed_database()
    return redirect(url_for('index'))

@app.route('/hub')
def hub(): return render_template('layout.html', mode='hub', links=RESEARCH_LINKS, app_name=APP_NAME)

@app.route('/tts', methods=['GET', 'POST'])
def tts_route():
    if request.method == 'POST':
        text = request.form.get('text', '')
        accent = request.form.get('accent', 'com')
    else:
        text = request.args.get('text', '')
        accent = request.args.get('accent', 'com')
    fp = BytesIO()
    try: 
        gTTS(text=text, lang='en', tld=accent).write_to_fp(fp)
        fp.seek(0)
        return send_file(fp, mimetype='audio/mpeg')
    except: return "Error", 500

@app.route('/systems', methods=['GET', 'POST'])
def systems():
    analysis = None
    if request.method == 'POST':
        w = request.form.get('word', '').strip().lower()
        if w:
            morph = [f"<b>{k}</b>: {v}" for k,v in utils.MORPHOLOGY_DB.items() if k in w]
            try: tr = GoogleTranslator(source='en', target='pt').translate(w)
            except: tr = "..."
            ipa_val = utils.get_phonetic(w) 
            analysis = {"word": w, "morph": morph, "trans": tr, "ipa": ipa_val}
    return render_template('layout.html', mode='systems', analysis=analysis, app_name=APP_NAME)

@app.route('/phrases')
def phrases(): 
    # Passamos o banco de frases aqui
    return render_template('layout.html', mode='phrases', phrases=ACADEMIC_PHRASEBANK, app_name=APP_NAME)

@app.route('/miner', methods=['GET', 'POST'])
def miner():
    if request.method == 'POST':
        full_text = request.form.get('texto_full', '')
        for f in request.files.getlist("arquivo_upload"):
            try:
                if f.filename.endswith('.pdf'): pdf = PdfReader(f); full_text += " ".join([page.extract_text() or "" for page in pdf.pages])
                elif f.filename.endswith('.txt'): full_text += f.read().decode('utf-8', errors='ignore')
                elif f.filename.endswith(('.ris','.nbib')): full_text += utils.parse_ris_nbib(f.read().decode('utf-8', errors='ignore'))
            except: pass
        return render_template('layout.html', mode='miner_res', keywords=utils.get_top_keywords(full_text), app_name=APP_NAME)
    return render_template('layout.html', mode='miner_home', app_name=APP_NAME)

@app.route('/search', methods=['GET', 'POST'])
def search():
    results = []
    query = ""
    if request.method == 'POST':
        query = request.form.get('query', '')
        if query:
            try:
                handle = Entrez.esearch(db="pubmed", term=query, retmax=5) 
                record = Entrez.read(handle)
                id_list = record["IdList"]
                if id_list:
                    handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
                    raw_data = handle.read()
                    articles = utils.parse_nbib_bulk(raw_data)
                    for art in articles:
                        if art.get('abstract'):
                            results.append({"title": art['title'],"abstract": art['abstract'],"journal": "PubMed Article"})
            except: pass
    return render_template('layout.html', mode='search', results=results, query=query, app_name=APP_NAME)

@app.route('/summarizer', methods=['GET', 'POST'])
def summarizer():
    batch_results = []
    if request.method == 'POST':
        content_to_process = []
        if 'nbib_file' in request.files:
            f = request.files['nbib_file']
            if f.filename:
                try:
                    content = f.read().decode('utf-8', errors='replace')
                    all_articles = utils.parse_nbib_bulk(content)
                    content_to_process = all_articles[:5] 
                except Exception as e: print(e)
        text_input = request.form.get('text_input', '')
        title_input = request.form.get('title_input', 'Manual Input')
        if text_input: content_to_process.append({"title": title_input, "abstract": text_input})
            
        for art in content_to_process:
            if art['abstract']:
                formatted_html = utils.format_abstract_smart(art['abstract'])
                clean_text_for_api = re.sub(r'<.*?>', ' ', formatted_html)
                try:
                    time.sleep(0.5)
                    translation = GoogleTranslator(source='en', target='pt').translate(clean_text_for_api[:4500])
                except: translation = "Tradução indisponível."
                batch_results.append({
                    "title": art['title'], 
                    "formatted_html": formatted_html, 
                    "clean_text": clean_text_for_api, 
                    "translation": translation
                })
        if batch_results: log_activity(len(batch_results) * 2)
    return render_template('layout.html', mode='summarizer', batch_results=batch_results, app_name=APP_NAME)

@app.route('/get_quiz', methods=['POST'])
def get_quiz():
    text_content = request.form.get('text', '')
    clean_text = re.sub(r'<[^>]*>', ' ', text_content)
    clean_text = re.sub(r'\s+', ' ', clean_text).strip()
    quiz_html = utils.generate_smart_quiz(clean_text)
    return jsonify({"html": quiz_html})

@app.route('/pronunciation')
def pronunciation():
    words = [c.front for c in Card.query.all()]
    if not words: words = ["Microbiology", "Respiration", "Analysis", "Data", "Science", "Health"]
    word_data = []
    seen = set()
    for w in words:
        clean = re.sub(r'[^\w\s]', '', w).strip()
        if len(clean) > 3 and clean.lower() not in seen:
            word_data.append({"w": clean, "ipa": utils.get_phonetic(clean)})
            seen.add(clean.lower())
    return render_template('layout.html', mode='pronunciation', words_json=json.dumps(word_data), app_name=APP_NAME)

@app.route('/vocab_list')
def vocab_list():
    cards = Card.query.order_by(Card.id.desc()).all()
    return render_template('layout.html', mode='vocab_list', cards=cards, app_name=APP_NAME)

@app.route('/export_vocab')
def export_vocab():
    cards = Card.query.all()
    si = StringIO()
    cw = csv.writer(si)
    cw.writerow(['Front', 'Back', 'IPA', 'Context', 'Next Review'])
    for c in cards:
        cw.writerow([c.front, c.back, c.ipa, c.context or "", c.next_review])
    output = BytesIO()
    output.write(u'\ufeff'.encode('utf-8'))
    output.write(si.getvalue().encode('utf-8'))
    output.seek(0)
    return send_file(output, mimetype='text/csv', as_attachment=True, download_name=f'scifluency_{datetime.now().strftime("%Y-%m-%d")}.csv')

@app.route('/import_vocab', methods=['POST'])
def import_vocab():
    if 'csv_file' in request.files:
        f = request.files['csv_file']
        if f.filename:
            try:
                stream = TextIOWrapper(f.stream, encoding='utf-8')
                reader = csv.reader(stream)
                header = next(reader, None)
                count = 0
                for row in reader:
                    if len(row) >= 2:
                        front = row[0]; back = row[1]
                        ipa_text = row[2] if len(row) > 2 else utils.get_phonetic(front)
                        context_text = row[3] if len(row) > 3 else ""
                        exists = Card.query.filter_by(front=front, deck_id="my_vocab").first()
                        if not exists:
                            db.session.add(Card(front=front, back=back, ipa=ipa_text, context=context_text, deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d')))
                            count += 1
                db.session.commit()
            except Exception as e: print(e)
    return redirect(url_for('vocab_list'))

@app.route('/delete_card/<int:id>')
def delete_card(id):
    c = Card.query.get(id)
    if c:
        db.session.delete(c); db.session.commit()
    return redirect(url_for('vocab_list'))

@app.route('/adicionar_vocab', methods=['POST'])
def adicionar_vocab():
    original = request.form.get('term')
    context_text = request.form.get('context', '')
    is_ajax = request.form.get('ajax', '0') == '1'
    
    if original:
        en, pt = utils.process_language_logic(original)
        front = en; back = pt
    else:
        front = request.form.get('f'); back = request.form.get('b')
    
    ipa_text = utils.get_phonetic(front)
    db.session.add(Card(front=front, back=back, ipa=ipa_text, context=context_text, deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d'), interval=0))
    db.session.commit()
    log_activity(1)
    
    if is_ajax: return jsonify({"status": "success", "msg": "Salvo!"})
    return redirect(request.referrer)

@app.route('/jogar/<id>')
def jogar(id):
    target_deck = "my_vocab"
    due = Card.query.filter(Card.deck_id == target_deck, Card.next_review <= datetime.now().strftime('%Y-%m-%d')).all()
    card_data = None
    if due:
        c = random.choice(due)
        cloze_parts = re.split(r'\[(.*?)\]', c.front) 
        is_cloze = len(cloze_parts) > 1
        card_data = {"id": c.id, "front": c.front, "back": c.back, "ipa": c.ipa, "context": c.context, "is_cloze": is_cloze, "cloze_hint": (cloze_parts[0] + "_____" + cloze_parts[2]) if is_cloze else ""}
    return render_template('layout.html', mode='study_play', deck=Deck.query.get(target_deck), card=card_data, app_name=APP_NAME)

@app.route('/rate_card', methods=['POST'])
def rate_card():
    c = Card.query.filter_by(deck_id="my_vocab", front=request.form.get('front')).first()
    if c:
        rating = request.form.get('rating')
        new_int = 1 if rating == 'hard' else (int((c.interval or 0) * 1.5) + 1 if rating == 'medium' else int((c.interval or 0) * 2.5) + 2)
        c.interval = new_int
        c.next_review = (datetime.now() + timedelta(days=new_int)).strftime('%Y-%m-%d')
        db.session.commit()
        log_activity(1)
    return redirect(url_for('jogar', id="my_vocab"))

@app.route('/tutor', methods=['GET', 'POST'])
def tutor():
    res = None
    if request.method == 'POST':
        user_input = request.form.get('pt_text', '')
        if user_input:
            en_text, pt_text = utils.process_language_logic(user_input)
            words = re.findall(r'\b[a-zA-Z]{4,}\b', en_text)
            vocab_breakdown = []
            unique_words = set([w.lower() for w in words])
            for w in unique_words:
                tr = GoogleTranslator(source='en', target='pt').translate(w)
                ipa_val = utils.get_phonetic(w) 
                vocab_breakdown.append({"word": w, "trans": tr, "ipa": ipa_val})
            res = {"pt": pt_text, "en": en_text, "vocab": vocab_breakdown}
    return render_template('layout.html', mode='tutor', res=res, app_name=APP_NAME)

@app.route('/checker', methods=['GET', 'POST'])
def checker():
    original = ""
    corrected = ""
    if request.method == 'POST':
        original = request.form.get('text_input', '')
        if original:
            corrected = utils.improve_english_text(original)
    return render_template('layout.html', mode='checker', original=original, corrected=corrected, app_name=APP_NAME)

@app.route('/processar', methods=['POST'])
def processar():
    text_content = request.form.get('texto_full', '')
    if 'arquivo_upload' in request.files:
        f = request.files['arquivo_upload']
        try:
            if f.filename.endswith('.pdf'): text_content += " ".join([p.extract_text() or "" for p in PdfReader(f).pages])
            elif f.filename.endswith('.txt'): text_content += f.read().decode('utf-8', errors='ignore')
            elif f.filename.endswith(('.ris','.nbib')): text_content += utils.parse_ris_nbib(f.read().decode('utf-8',errors='ignore'))
        except: pass
    clean_text = re.sub(r'\s+', ' ', text_content)
    frases = re.split(r'(?<=[.!?])\s+(?=[A-Z])', clean_text)
    sid = str(uuid.uuid4())[:8]
    db.session.add(Story(id=sid, title=(frases[0][:40] + "..." if frases else "New Text")))
    tr = GoogleTranslator(source='en', target='pt')
    for f in frases[:50]:
        if len(f)<15: continue
        try: pt=tr.translate(f)
        except: pt="..."
        db.session.add(Sentence(en=f, pt=pt, story_id=sid))
    db.session.commit()
    log_activity(5)
    return redirect(url_for('ler', id=sid))

@app.route('/ler/<id>')
def ler(id): return render_template('layout.html', mode='read', story=Story.query.get(id), app_name=APP_NAME)

@app.route('/deletar/<id>')
def deletar(id): 
    db.session.delete(Story.query.get(id)); db.session.commit()
    return redirect(url_for('index'))

@app.route('/traduzir_palavra')
def traduzir_palavra(): return jsonify({"t": GoogleTranslator(source='en', target='pt').translate(request.args.get('w',''))})

@app.route('/novo')
def novo(): return render_template('layout.html', mode='new', app_name=APP_NAME)

@app.route('/podcast/<id>')
def podcast(id):
    story = Story.query.get(id)
    if not story: return "Texto não encontrado", 404
    
    try:
        full_audio = BytesIO()
        intro_text = f"SciFluency Bilingual Audio. {story.title}."
        tts_intro = gTTS(text=intro_text, lang='en', tld='com')
        tts_intro.write_to_fp(full_audio)
        
        count = 0
        for s in story.sentences:
            if count > 20: break 
            if s.en and s.pt:
                tts_en = gTTS(text=s.en, lang='en', tld='com')
                tts_en.write_to_fp(full_audio)
                tts_pt = gTTS(text=s.pt, lang='pt', tld='com.br')
                tts_pt.write_to_fp(full_audio)
                time.sleep(0.3)
            count += 1
            
        tts_end = gTTS(text="End of session.", lang='en')
        tts_end.write_to_fp(full_audio)
        full_audio.seek(0)
        return send_file(full_audio, mimetype='audio/mpeg', as_attachment=True, download_name=f'bilingual_{id}.mp3')

    except Exception as e:
        print(f"Erro Podcast: {e}")
        return f"Erro ao gerar podcast: {e}", 500

with app.app_context(): 
    db.create_all()
    check_and_migrate_db()
    seed_database()

if __name__ == '__main__': app.run(host='0.0.0.0', port=5000)
