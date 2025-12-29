import os
import random
import re
import time
import json
import csv
from io import BytesIO, StringIO, TextIOWrapper
from datetime import datetime, timedelta
from collections import Counter

from flask import Flask, render_template_string, request, redirect, url_for, send_file, jsonify
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import text, inspect
from deep_translator import GoogleTranslator
from pypdf import PdfReader
from gtts import gTTS
from Bio import Entrez
import eng_to_ipa as ipa

# --- CONFIGURAÇÃO ---
app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = os.environ.get('DATABASE_URL', 'sqlite:///scifluency.db')
if app.config['SQLALCHEMY_DATABASE_URI'].startswith("postgres://"):
    app.config['SQLALCHEMY_DATABASE_URI'] = app.config['SQLALCHEMY_DATABASE_URI'].replace("postgres://", "postgresql://", 1)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

Entrez.email = "researcher@example.com"
Entrez.tool = "SciFluencyResearch"

db = SQLAlchemy(app)
APP_NAME = "SciFluency"

# --- MODELOS ---
class Deck(db.Model):
    id = db.Column(db.String(50), primary_key=True)
    name = db.Column(db.String(100))
    pt_name = db.Column(db.String(100))
    icon = db.Column(db.String(10))
    cards = db.relationship('Card', backref='deck', lazy=True, cascade="all, delete-orphan")

class Card(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    front = db.Column(db.String(200)) 
    back = db.Column(db.String(200))
    ipa = db.Column(db.String(200))
    context = db.Column(db.Text)
    deck_id = db.Column(db.String(50), db.ForeignKey('deck.id'))
    next_review = db.Column(db.String(20))
    interval = db.Column(db.Integer, default=0)

class Story(db.Model):
    id = db.Column(db.String(50), primary_key=True)
    title = db.Column(db.String(200))
    sentences = db.relationship('Sentence', backref='story', lazy=True, cascade="all, delete-orphan")

class Sentence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    en = db.Column(db.Text)
    pt = db.Column(db.Text)
    story_id = db.Column(db.String(50), db.ForeignKey('story.id'))

class StudyLog(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    date = db.Column(db.String(10))
    count = db.Column(db.Integer, default=0)

# --- DADOS E CONSTANTES ---
RESEARCH_LINKS = [
    {"name": "PubMed", "url": "https://pubmed.ncbi.nlm.nih.gov/", "icon": "🧬", "desc": "Biomedical Literature"},
    {"name": "SciELO", "url": "https://scielo.org/", "icon": "🌎", "desc": "Open Access Journals"},
    {"name": "Google Scholar", "url": "https://scholar.google.com.br/", "icon": "🎓", "desc": "Academic Search"},
    {"name": "Cochrane", "url": "https://www.cochranelibrary.com/", "icon": "🏥", "desc": "Evidence-Based Medicine"},
    {"name": "Scopus", "url": "https://www.scopus.com/", "icon": "🔭", "desc": "Citation Database"},
    {"name": "Web of Science", "url": "https://www.webofscience.com/", "icon": "🕸️", "desc": "Scientific Citation Index"}
]

ACADEMIC_PHRASEBANK = {
    "1. Introduction & Context": [
        {"en": "Recent developments in this field have heightened the need for...", "pt": "Desenvolvimentos recentes neste campo aumentaram a necessidade de..."},
        {"en": "Currently, there is a paucity of data regarding...", "pt": "Atualmente, há escassez de dados sobre..."},
        {"en": "This study aims to investigate the relationship between...", "pt": "Este estudo visa investigar a relação entre..."},
        {"en": "Previous research has established that...", "pt": "Pesquisas anteriores estabeleceram que..."},
        {"en": "The primary objective of this paper is to evaluate...", "pt": "O objetivo principal deste artigo é avaliar..."}
    ],
    "2. Methods & Materials": [
        {"en": "Data were collected using a semi-structured interview guide.", "pt": "Os dados foram coletados usando um roteiro de entrevista semiestruturado."},
        {"en": "The participants were divided into two groups.", "pt": "Os participantes foram divididos em dois grupos."},
        {"en": "Statistical analysis was performed using SPSS software.", "pt": "A análise estatística foi realizada usando o software SPSS."},
        {"en": "We excluded patients with a history of...", "pt": "Excluímos pacientes com histórico de..."}
    ],
    "3. Results & Findings": [
        {"en": "There was a significant correlation between...", "pt": "Houve uma correlação significativa entre..."},
        {"en": "Table 1 presents the demographic characteristics of the sample.", "pt": "A Tabela 1 apresenta as características demográficas da amostra."},
        {"en": "The results indicate that...", "pt": "Os resultados indicam que..."},
        {"en": "Our findings are consistent with those of previous studies.", "pt": "Nossos achados são consistentes com os de estudos anteriores."}
    ],
    "4. Discussion & Argumentation": [
        {"en": "These findings suggest that...", "pt": "Esses achados sugerem que..."},
        {"en": "One possible explanation for this result is...", "pt": "Uma possível explicação para este resultado é..."},
        {"en": "This study provides new insights into...", "pt": "Este estudo fornece novos insights sobre..."},
        {"en": "However, some limitations should be noted.", "pt": "No entanto, algumas limitações devem ser notadas."}
    ],
    "5. Conclusion": [
        {"en": "In conclusion, this study demonstrates that...", "pt": "Em conclusão, este estudo demonstra que..."},
        {"en": "The evidence from this study suggests...", "pt": "As evidências deste estudo sugerem..."},
        {"en": "Overall, these results highlight the importance of...", "pt": "No geral, esses resultados destacam a importância de..."}
    ]
}

ACADEMIC_REPLACEMENTS = {
    "big": "substantial", "huge": "significant", "bad": "detrimental",
    "good": "beneficial", "think": "hypothesize", "get": "obtain",
    "make": "generate", "show": "demonstrate", "use": "utilize",
    "really": "significantly", "very": "highly", "look at": "examine",
    "prove": "validate", "change": "alter", "stop": "cease"
}
MANUAL_DICT = {
    "rins": "Kidneys", "rim": "Kidney", "pulmão": "Lung", "pulmoes": "Lungs",
    "coração": "Heart", "figado": "Liver", "cerebro": "Brain",
    "metástase": "Metastasis", "metastases": "Metastases", "casa": "House"
}
MORPHOLOGY_DB = {"un": "Not (Não)", "re": "Again (Novamente)", "itis": "Inflammation"}

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

# --- HELPERS ---
def get_phonetic(text):
    if not text: return ""
    words = text.split()
    ipa_result = []
    for word in words:
        clean = re.sub(r'[^\w\s]', '', word)
        if not clean: 
            ipa_result.append(word); continue
        try:
            generated = ipa.convert(clean)
            if "*" in generated: ipa_result.append(clean)
            else: ipa_result.append(generated)
        except: ipa_result.append(clean)
    return "/" + " ".join(ipa_result) + "/"

def process_language_logic(text):
    if not text: return None, None
    clean_text = text.strip()
    lower_text = clean_text.lower()
    if lower_text in MANUAL_DICT: return MANUAL_DICT[lower_text], clean_text
    try:
        trans_en = GoogleTranslator(source='pt', target='en').translate(clean_text)
        if trans_en and trans_en.lower() != lower_text: return trans_en, clean_text
        trans_pt = GoogleTranslator(source='en', target='pt').translate(clean_text)
        return clean_text, trans_pt
    except: return clean_text, "Error"

def improve_english_text(text):
    if not text: return ""
    try:
        pt = GoogleTranslator(source='en', target='pt').translate(text)
        fixed = GoogleTranslator(source='pt', target='en').translate(pt)
    except: fixed = text
    for s, a in ACADEMIC_REPLACEMENTS.items():
        fixed = re.sub(rf"\b{s}\b", a, fixed, flags=re.IGNORECASE)
    return fixed

def generate_smart_quiz(text):
    if not text: return ""
    words = text.split()
    output = []
    ignore = {"the","and","of","to","in","is","that","for","with","was","as","are","this","from","by"}
    candidates = []
    for i, w in enumerate(words):
        clean = re.sub(r'[^\w]', '', w).lower()
        if len(clean) > 4 and clean not in ignore:
            candidates.append(i)
    if candidates:
        to_hide = set(random.sample(candidates, k=max(1, int(len(candidates) * 0.25))))
    else:
        to_hide = set()
    for i, w in enumerate(words):
        if i in to_hide:
            output.append(f"<span class='quiz-hidden' onclick='reveal(this, \"{w}\")'>[ ? ]</span>")
        else:
            output.append(w)
    return " ".join(output)

def seed_database():
    if Deck.query.filter_by(id="my_vocab").first() is None:
        db.session.add(Deck(id="my_vocab", name="My Vocabulary", pt_name="Meu Banco de Palavras", icon="🗂️"))
        db.session.add(Card(front="Welcome", back="Bem-vindo", ipa="/ˈwel.kəm/", context="Welcome to SciFluency.", deck_id="my_vocab", next_review=datetime.now().strftime('%Y-%m-%d')))
        db.session.commit()

def parse_ris_nbib(content):
    text = ""
    for line in content.split('\n'):
        line = line.strip()
        if line.startswith(('TI', 'T1', 'AB')): parts = line.split('-', 1); 
        if len(parts) > 1: text += parts[1].strip() + "\n\n"
    return text

def get_top_keywords(text):
    words = re.findall(r'\b[a-zA-Z]{5,}\b', text.lower())
    stop = {"the","and","of","to","in","is","that","for","with","was","as","are","this","from","which", "study", "results", "using", "between", "during"}
    return Counter([w for w in words if w not in stop]).most_common(80)

def log_activity(points):
    today = datetime.now().strftime('%Y-%m-%d')
    log = StudyLog.query.filter_by(date=today).first()
    if not log: log = StudyLog(date=today, count=0); db.session.add(log)
    log.count += points
    db.session.commit()

def parse_nbib_bulk(content):
    articles = []
    current = {"title": "No Title", "abstract": ""}
    lines = content.replace('\r', '').split('\n')
    last_tag = None
    for line in lines:
        stripped_line = line.strip()
        if not stripped_line: continue
        if line.startswith("PMID-") or line.startswith("PMID -"):
            if current["abstract"] and len(current["abstract"]) > 50: articles.append(current)
            current = {"title": "No Title", "abstract": ""}
            last_tag = "PMID"
        elif line.startswith("TI  - "):
            current["title"] = line[6:].strip()
            last_tag = "TI"
        elif line.startswith("AB  - "):
            current["abstract"] = line[6:].strip()
            last_tag = "AB"
        elif line.startswith("      ") and last_tag == "AB":
            current["abstract"] += " " + line.strip()
        elif line.startswith("      ") and last_tag == "TI":
            current["title"] += " " + line.strip()
        elif re.match(r'^[A-Z]{2,4}\s*-', line):
            last_tag = None
    if current["abstract"] and len(current["abstract"]) > 50: articles.append(current)
    return articles

def format_abstract_smart(text):
    if not text: return ""
    formatted = text
    # Destaque de seções
    sections = ["BACKGROUND", "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSIONS", "CONCLUSION", "DISCUSSION"]
    for sec in sections:
        pattern = re.compile(rf"({sec}[:\s])", re.IGNORECASE)
        if pattern.search(formatted):
            formatted = pattern.sub(r"<br><br><b style='color:#2c3e50;'>\1</b>", formatted)
    
    sentences = re.split(r'(?<=[.!?])\s+', formatted)
    final_html = ""
    for sent in sentences:
        if len(sent.strip()) > 1:
            words_html = ""
            for word in sent.split():
                clean_w = re.sub(r"[^\w]", "", word) 
                if clean_w:
                    words_html += f"<span class='word-span' onclick='mineWord(event, \"{clean_w}\")'>{word}</span> "
                else:
                    words_html += word + " "
            # Usamos prepare(this) sem argumentos. O JS pega o texto do elemento.
            final_html += f"<div class='sentence-block' onclick='prepare(this)' style='margin-bottom:5px; padding:5px; cursor:pointer;'>{words_html}</div>"
    return final_html

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
    return render_template_string(PAGE_LAYOUT, mode='list', stories=Story.query.all(), decks=[deck_data], stats=stats, app_name=APP_NAME, links=RESEARCH_LINKS)

@app.route('/reset_decks')
def reset_decks():
    db.drop_all(); db.create_all(); seed_database()
    return redirect(url_for('index'))

@app.route('/hub')
def hub(): return render_template_string(PAGE_LAYOUT, mode='hub', links=RESEARCH_LINKS, app_name=APP_NAME)

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
        # OBTEM AUDIO DO GOOGLE (Servidor)
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
            morph = [f"<b>{k}</b>: {v}" for k,v in MORPHOLOGY_DB.items() if k in w]
            try: tr = GoogleTranslator(source='en', target='pt').translate(w)
            except: tr = "..."
            ipa_val = get_phonetic(w) 
            analysis = {"word": w, "morph": morph, "trans": tr, "ipa": ipa_val}
    return render_template_string(PAGE_LAYOUT, mode='systems', analysis=analysis, app_name=APP_NAME)

@app.route('/phrases')
def phrases(): return render_template_string(PAGE_LAYOUT, mode='phrases', phrases=ACADEMIC_PHRASEBANK, app_name=APP_NAME)

@app.route('/miner', methods=['GET', 'POST'])
def miner():
    if request.method == 'POST':
        full_text = request.form.get('texto_full', '')
        for f in request.files.getlist("arquivo_upload"):
            try:
                if f.filename.endswith('.pdf'): pdf = PdfReader(f); full_text += " ".join([page.extract_text() or "" for page in pdf.pages])
                elif f.filename.endswith('.txt'): full_text += f.read().decode('utf-8', errors='ignore')
                elif f.filename.endswith(('.ris','.nbib')): full_text += parse_ris_nbib(f.read().decode('utf-8', errors='ignore'))
            except: pass
        return render_template_string(PAGE_LAYOUT, mode='miner_res', keywords=get_top_keywords(full_text), app_name=APP_NAME)
    return render_template_string(PAGE_LAYOUT, mode='miner_home', app_name=APP_NAME)

@app.route('/search', methods=['GET', 'POST'])
def search():
    results = []
    query = ""
    retstart = 0
    count = 0
    
    if request.method == 'POST':
        query = request.form.get('query', '')
        try:
            retstart = int(request.form.get('retstart', 0))
        except:
            retstart = 0
            
        if query:
            try:
                handle = Entrez.esearch(db="pubmed", term=query, retmax=5, retstart=retstart) 
                record = Entrez.read(handle)
                count = int(record['Count'])
                id_list = record["IdList"]
                if id_list:
                    handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
                    raw_data = handle.read()
                    articles = parse_nbib_bulk(raw_data)
                    for art in articles:
                        if art.get('abstract'):
                            results.append({"title": art['title'],"abstract": art['abstract'],"journal": "PubMed Article"})
            except: pass
    return render_template_string(PAGE_LAYOUT, mode='search', results=results, query=query, retstart=retstart, total_count=count, app_name=APP_NAME)

@app.route('/save_article', methods=['POST'])
def save_article():
    title = request.form.get('title', 'No Title')
    abstract = request.form.get('abstract', '')
    if not abstract:
        return redirect(url_for('search'))
    sid = str(uuid.uuid4())[:8]
    new_story = Story(id=sid, title=title[:150]) 
    db.session.add(new_story)
    clean_text = re.sub(r'\s+', ' ', abstract)
    frases = re.split(r'(?<=[.!?])\s+(?=[A-Z])', clean_text)
    tr = GoogleTranslator(source='en', target='pt')
    for f in frases:
        if len(f) < 10: continue
        try: pt = tr.translate(f)
        except: pt = "..."
        db.session.add(Sentence(en=f, pt=pt, story_id=sid))
    db.session.commit()
    log_activity(3)
    return redirect(url_for('index'))

@app.route('/summarizer', methods=['GET', 'POST'])
def summarizer():
    batch_results = []
    if request.method == 'POST':
        content_to_process = []
        if 'nbib_file' in request.files:
            f = request.files['nbib_file']
            if f.filename:
                try:
                    if f.filename.lower().endswith('.pdf'):
                        raw_pdf_text = " ".join([page.extract_text() or "" for page in PdfReader(f).pages])
                        clean_pdf_text = re.sub(r'\s+', ' ', raw_pdf_text).strip()
                        content_to_process = [{"title": f.filename, "abstract": clean_pdf_text}]
                    else:
                        content = f.read().decode('utf-8', errors='replace')
                        all_articles = parse_nbib_bulk(content)
                        content_to_process = all_articles[:5] 
                except Exception as e: 
                    print(f"Erro no Summarizer: {e}")

        text_input = request.form.get('text_input', '')
        title_input = request.form.get('title_input', 'Manual Input')
        if text_input: content_to_process.append({"title": title_input, "abstract": text_input})
            
        for art in content_to_process:
            if art['abstract']:
                formatted_html = format_abstract_smart(art['abstract'])
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
    return render_template_string(PAGE_LAYOUT, mode='summarizer', batch_results=batch_results, app_name=APP_NAME)

@app.route('/get_quiz', methods=['POST'])
def get_quiz():
    text_content = request.form.get('text', '')
    clean_text = re.sub(r'<[^>]*>', ' ', text_content)
    clean_text = re.sub(r'\s+', ' ', clean_text).strip()
    quiz_html = generate_smart_quiz(clean_text)
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
            word_data.append({"w": clean, "ipa": get_phonetic(clean)})
            seen.add(clean.lower())
    return render_template_string(PAGE_LAYOUT, mode='pronunciation', words_json=json.dumps(word_data), app_name=APP_NAME)

@app.route('/vocab_list')
def vocab_list():
    cards = Card.query.order_by(Card.id.desc()).all()
    return render_template_string(PAGE_LAYOUT, mode='vocab_list', cards=cards, app_name=APP_NAME)

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
                        ipa_text = row[2] if len(row) > 2 else get_phonetic(front)
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
        db.session.delete(c)
        db.session.commit()
    return redirect(url_for('vocab_list'))

@app.route('/adicionar_vocab', methods=['POST'])
def adicionar_vocab():
    original = request.form.get('term')
    context_text = request.form.get('context', '')
    is_ajax = request.form.get('ajax', '0') == '1'
    if original:
        en, pt = process_language_logic(original)
        front = en
        back = pt
    else:
        front = request.form.get('f')
        back = request.form.get('b')
    ipa_text = get_phonetic(front)
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
    return render_template_string(PAGE_LAYOUT, mode='study_play', deck=Deck.query.get(target_deck), card=card_data, app_name=APP_NAME)

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
            en_text, pt_text = process_language_logic(user_input)
            words = re.findall(r'\b[a-zA-Z]{4,}\b', en_text)
            vocab_breakdown = []
            unique_words = set([w.lower() for w in words])
            for w in unique_words:
                tr = GoogleTranslator(source='en', target='pt').translate(w)
                ipa_val = get_phonetic(w) 
                vocab_breakdown.append({"word": w, "trans": tr, "ipa": ipa_val})
            res = {"pt": pt_text, "en": en_text, "vocab": vocab_breakdown}
    return render_template_string(PAGE_LAYOUT, mode='tutor', res=res, app_name=APP_NAME)

@app.route('/checker', methods=['GET', 'POST'])
def checker():
    original = ""
    corrected = ""
    if request.method == 'POST':
        original = request.form.get('text_input', '')
        if original:
            corrected = improve_english_text(original)
    return render_template_string(PAGE_LAYOUT, mode='checker', original=original, corrected=corrected, app_name=APP_NAME)

@app.route('/processar', methods=['POST'])
def processar():
    text_content = request.form.get('texto_full', '')
    uploaded_files = request.files.getlist("arquivo_upload")
    processed_count = 0
    last_story_id = None

    if text_content.strip():
        sid = str(uuid.uuid4())[:8]
        clean_text = re.sub(r'\s+', ' ', text_content)
        frases = re.split(r'(?<=[.!?])\s+(?=[A-Z])', clean_text)
        db.session.add(Story(id=sid, title=(frases[0][:40] + "..." if frases else "Pasted Text")))
        tr = GoogleTranslator(source='en', target='pt')
        for f in frases[:50]:
            if len(f)<10: continue
            try: pt=tr.translate(f)
            except: pt="..."
            db.session.add(Sentence(en=f, pt=pt, story_id=sid))
        processed_count += 1
        last_story_id = sid

    for f in uploaded_files:
        if not f or not f.filename: continue
        file_text = ""
        try:
            if f.filename.endswith('.pdf'): 
                file_text = " ".join([p.extract_text() or "" for p in PdfReader(f).pages])
            elif f.filename.endswith('.txt'): 
                file_text = f.read().decode('utf-8', errors='ignore')
            elif f.filename.endswith(('.ris','.nbib')): 
                file_text = parse_ris_nbib(f.read().decode('utf-8',errors='ignore'))
        except: continue

        if file_text.strip():
            clean_text = re.sub(r'\s+', ' ', file_text)
            frases = re.split(r'(?<=[.!?])\s+(?=[A-Z])', clean_text)
            sid = str(uuid.uuid4())[:8]
            title = f.filename
            if len(frases) > 0: title = f.filename + " - " + frases[0][:30]
            db.session.add(Story(id=sid, title=title))
            tr = GoogleTranslator(source='en', target='pt')
            for sentence in frases[:60]:
                if len(sentence) < 15: continue
                try: pt = tr.translate(sentence)
                except: pt = "..."
                db.session.add(Sentence(en=sentence, pt=pt, story_id=sid))
            processed_count += 1
            last_story_id = sid

    db.session.commit()
    log_activity(5 * processed_count)
    if processed_count == 1 and last_story_id:
        return redirect(url_for('ler', id=last_story_id))
    return redirect(url_for('index'))

@app.route('/ler/<id>')
def ler(id): return render_template_string(PAGE_LAYOUT, mode='read', story=Story.query.get(id), app_name=APP_NAME)

@app.route('/deletar/<id>')
def deletar(id): 
    db.session.delete(Story.query.get(id)); db.session.commit()
    return redirect(url_for('index'))

@app.route('/traduzir_palavra')
def traduzir_palavra(): return jsonify({"t": GoogleTranslator(source='en', target='pt').translate(request.args.get('w',''))})

@app.route('/novo')
def novo(): return render_template_string(PAGE_LAYOUT, mode='new', app_name=APP_NAME)

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

# --- FRONTEND (HÍBRIDO: AUDIO SERVER + JS TIMER) ---
PAGE_LAYOUT = r"""
<!DOCTYPE html><html><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{{ app_name }}</title>
<style>
:root { --bg: #f8faff; --card: #fff; --text: #333; --subtext: #7f8c8d; --accent: #0088cc; --border: #e1e8ed; --input: #fff; --heat-0: #ecf0f1; --heat-1: #a8e6cf; --heat-2: #55efc4; --heat-3: #00b894; --shadow: rgba(0,0,0,0.03); }
body.dark-mode { --bg: #121212; --card: #1e1e1e; --text: #e0e0e0; --subtext: #b0b0b0; --accent: #4fc3f7; --border: #333; --input: #2c2c2c; --heat-0: #2c2c2c; --heat-1: #1e3a2f; --heat-2: #2e7d62; --heat-3: #00b894; --shadow: rgba(0,0,0,0.5); }
body{font-family:'Segoe UI',sans-serif;background:var(--bg);margin:0;padding-bottom:180px;color:var(--text); transition: background 0.3s;}
.header{background:var(--card);padding:18px;border-bottom:1px solid var(--border);display:flex;align-items:center;position:sticky;top:0;z-index:1000; justify-content: space-between;}
.container{padding:20px;max-width:650px;margin:0 auto;}
.sidebar{height:100%;width:0;position:fixed;top:0;left:0;background:var(--card);overflow-x:hidden;transition:0.3s;padding-top:60px;z-index:2000;box-shadow:2px 0 100px rgba(0,0,0,0.2);}
.sidebar a{padding:15px 25px;text-decoration:none;font-size:1.1rem;color:var(--text);display:block;border-bottom:1px solid var(--border);}
.sub-text {display:block;font-size:0.85rem;color:var(--subtext);font-weight:normal;margin-top:3px;}
.sidebar a span{display:block;font-size:0.8rem;color:var(--subtext);margin-top:2px;}
.card{background:var(--card);padding:20px;margin-bottom:15px;border-radius:16px;border:1px solid var(--border);box-shadow:0 4px 10px var(--shadow);}
.btn{background:var(--accent);color:#fff;border:none;padding:12px 20px;border-radius:25px;font-weight:bold;cursor:pointer;}
.sentence-block{padding:18px;margin-bottom:12px;border-left:5px solid transparent;cursor:pointer;background:var(--card);border-radius:8px; position: relative;}
.k-sent { cursor: pointer; transition: background-color 0.2s; padding: 2px 0; border-radius: 4px; }
.k-sent:hover { background-color: #e1f5fe; color: #000; }
.word-active { background-color: #ffeb3b !important; color: #000; box-shadow: 0 0 5px #ffeb3b; }
.player{position:fixed;bottom:0;left:0;width:100%;background:var(--card);border-top:1px solid var(--border);padding:20px;z-index:1500;display:flex;flex-direction:column;gap:10px;align-items:center;}
.controls-row{display:flex;gap:15px;align-items:center;}
input,textarea,select{width:100%;padding:12px;margin-bottom:10px;border:1px solid var(--border);border-radius:10px; background:var(--input); color:var(--text);}
.sys-box{background:var(--bg);padding:10px;border-radius:8px;margin-bottom:10px;}
.hub-card{text-align:center;padding:20px;cursor:pointer;}.hub-card:hover{filter:brightness(0.95);}
.badge {background:#ff7675;color:white;padding:2px 8px;border-radius:12px;font-size:0.8rem;margin-left:5px;}
.srs-btn {flex:1; border:none; padding:15px; border-radius:12px; font-weight:bold; cursor:pointer; color:#444;}
.modal-actions { display: flex; gap: 5px; }
.checker-box { display: flex; flex-direction: column; gap: 10px; }
.checker-label { font-weight: bold; color: var(--subtext); font-size: 0.9rem; margin-bottom: 5px; display: block;}
.text-area-box { width: 100%; padding: 15px; border-radius: 12px; border: 1px solid var(--border); font-family: inherit; font-size: 1rem; line-height: 1.5; resize: vertical; min-height: 120px; box-sizing: border-box; background:var(--input); color:var(--text); }
.corrected-box { border: 1px solid #2ecc71; background: #f0fdf4; color: #14532d; }
.arrow-separator { text-align: center; font-size: 1.5rem; color: var(--accent); margin: 5px 0; }
.save-form { margin-top: 15px; padding: 15px; background: var(--bg); border: 1px dashed #2ecc71; border-radius: 10px; }
.fab { position: fixed; bottom: 90px; right: 20px; width: 56px; height: 56px; background: var(--accent); color: white; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 28px; box-shadow: 0 4px 10px rgba(0,0,0,0.3); cursor: pointer; z-index: 1400; transition: transform 0.2s; user-select: none; }
.fab:active { transform: scale(0.95); }
.mic-btn-active { background-color: #e74c3c !important; animation: pulse 1.5s infinite; }
@keyframes pulse { 0% { box-shadow: 0 0 0 0 rgba(231, 76, 60, 0.7); } 70% { box-shadow: 0 0 0 10px rgba(231, 76, 60, 0); } 100% { box-shadow: 0 0 0 0 rgba(231, 76, 60, 0); } }
.dash-row { display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin-bottom: 15px; }
.stat-card { background: var(--card); padding: 15px; border-radius: 12px; text-align: center; box-shadow: 0 2px 5px var(--shadow); }
.stat-val { font-size: 1.8rem; font-weight: bold; color: var(--accent); }
.stat-lbl { font-size: 0.8rem; color: var(--subtext); text-transform: uppercase; letter-spacing: 1px; }
.heatmap { display: grid; grid-template-columns: repeat(7, 1fr); gap: 5px; margin-top: 10px; }
.heat-day { aspect-ratio: 1; border-radius: 4px; display: flex; align-items: center; justify-content: center; font-size: 0.7rem; color: #555; }
.cloze-blank { border-bottom: 2px solid var(--accent); color: transparent; background: #e3f2fd; padding: 0 5px; border-radius: 3px; }
.cloze-reveal .cloze-blank { color: var(--accent); background: transparent; font-weight: bold; }
.modal { display: none; position: fixed; z-index: 3000; left: 0; top: 0; width: 100%; height: 100%; overflow: auto; background-color: rgba(0,0,0,0.6); align-items: center; justify-content: center; }
.modal form { background-color: var(--card); margin: 15% auto; padding: 25px; border: 1px solid var(--border); width: 85%; max-width: 400px; border-radius: 16px; box-shadow: 0 10px 30px rgba(0,0,0,0.3); color: var(--text); position: relative; }
.word-row { display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid var(--border); padding:10px 0; }
.shadow-btn { float: right; margin-left: 10px; width: 40px; height: 40px; border-radius: 50%; padding: 0; display: flex; align-items: center; justify-content: center; font-size: 1.2rem; background: #ecf0f1; color: #333; }
.shadow-recording { background: #e74c3c !important; color: white; animation: pulse 1.5s infinite; }

.word-span { transition: 0.1s; border-bottom: 1px dotted transparent; }
body.mode-read .word-span { pointer-events: none; } 
body.mode-read .k-sent { cursor: pointer; }
body.mode-mine .k-sent { pointer-events: none; }
body.mode-mine .word-span { pointer-events: auto; cursor: context-menu; border-bottom: 2px dotted #e67e22; color: #d35400; }
body.mode-mine .word-span:hover { background: #f39c12; color: white; border-radius: 3px; }

/* QUIZ STYLES */
.quiz-hidden { background: #f1c40f; color: transparent; border-radius: 4px; padding: 0 5px; cursor: pointer; user-select: none; border-bottom: 2px solid #d35400; font-size: 0.9em; min-width: 30px; display: inline-block; text-align: center; }
.quiz-hidden:hover { background: #f39c12; }
.quiz-revealed { background: transparent; color: #2ecc71; font-weight: bold; border-bottom: none; animation: pop 0.3s ease; }
@keyframes pop { 0% { transform: scale(0.8); } 50% { transform: scale(1.1); } 100% { transform: scale(1); } }
</style>
<script>
let curAud = null;
let currentSpans = [];
let karaokeInterval = null;
let rate = 1.0;
let acc = 'com'; 
let shadowRec = null; 
let shadowChunks = [];

function openNav(){document.getElementById("side").style.width="280px";}
function closeNav(){document.getElementById("side").style.width="0";}
function toggleTheme() { document.body.classList.toggle('dark-mode'); const isDark = document.body.classList.contains('dark-mode'); localStorage.setItem('theme', isDark ? 'dark' : 'light'); document.getElementById('themeIcon').innerText = isDark ? '☀️' : '🌙'; }
window.onload = () => { document.body.classList.add('mode-read'); if(localStorage.getItem('theme') === 'dark') { document.body.classList.add('dark-mode'); if(document.getElementById('themeIcon')) document.getElementById('themeIcon').innerText = '☀️'; } };

// --- ESTRATÉGIA HÍBRIDA: ÁUDIO DO SERVIDOR + KARAOKE MATEMÁTICO ---
async function speak(text){ 
    if(curAud) { curAud.pause(); curAud = null; } 
    if(karaokeInterval) clearInterval(karaokeInterval);
    document.querySelectorAll('.word-active').forEach(w=>w.classList.remove('word-active'));

    const btn = document.getElementById('playBtn'); btn.innerText='⏳'; 
    const cleanTxt = text.replace(/<[^>]*>/g, "").replace(/\[|\]/g, "");

    try { 
        const response = await fetch('/tts', { 
            method: 'POST', 
            headers: { 'Content-Type': 'application/x-www-form-urlencoded' }, 
            body: 'text=' + encodeURIComponent(cleanTxt) + '&accent=' + acc 
        }); 
        
        const blob = await response.blob(); 
        const url = URL.createObjectURL(blob); 
        curAud = new Audio(url); 
        curAud.playbackRate = rate; 
        
        curAud.onloadedmetadata = () => { 
            btn.innerText='⏸';
            curAud.play();
            // Inicia o Karaoke Estimado
            startEstimatedKaraoke(curAud.duration);
        }; 
        
        curAud.onended = () => { 
            btn.innerText='▶'; 
            if(karaokeInterval) clearInterval(karaokeInterval);
            document.querySelectorAll('.word-active').forEach(w=>w.classList.remove('word-active')); 
        }; 
        
    } catch (error) { console.error("TTS Error:", error); btn.innerText='⚠️'; } 
}

function startEstimatedKaraoke(duration) {
    if (currentSpans.length === 0) return;
    
    // Tempo total em ms dividido pelo número de palavras
    // Ajustado pela velocidade de reprodução (rate)
    let totalTimeMs = (duration * 1000) / rate;
    let timePerWord = totalTimeMs / currentSpans.length;
    
    let index = 0;
    
    // Função para acender a próxima palavra
    const tick = () => {
        if (index >= currentSpans.length) {
            clearInterval(karaokeInterval);
            return;
        }
        
        // Remove destaque anterior
        if (index > 0) currentSpans[index-1].classList.remove('word-active');
        
        // Adiciona destaque atual
        currentSpans[index].classList.add('word-active');
        index++;
    };
    
    // Executa imediatamente e depois a cada intervalo
    tick();
    karaokeInterval = setInterval(tick, timePerWord);
}

function togglePlay(){ 
    if(!curAud) return; 
    const btn=document.getElementById('playBtn'); 
    if(curAud.paused){ 
        curAud.play(); 
        btn.innerText='⏸';
        // Recalcular karaoke seria complexo aqui, então simplificamos:
        // Se pausar, o karaoke pode dessincronizar levemente nesta versão simples
    } else{ 
        curAud.pause(); 
        btn.innerText='▶'; 
        if(karaokeInterval) clearInterval(karaokeInterval);
    } 
}

function prepare(el){ 
    if(curAud){curAud.pause(); if(karaokeInterval) clearInterval(karaokeInterval);} 
    document.querySelectorAll('.word-active').forEach(w=>w.classList.remove('word-active')); 
    
    // Tenta pegar spans de palavras individuais
    currentSpans = Array.from(el.querySelectorAll('.word-span')); 
    
    // Se não achou spans (ex: modo leitura simples), tenta spans de frases
    if (currentSpans.length === 0) currentSpans = Array.from(el.querySelectorAll('.k-sent'));
    
    // Se ainda vazio, usa o próprio elemento clicado como bloco único
    if (currentSpans.length === 0) currentSpans = [el];
    
    let textToRead = el.innerText; 
    speak(textToRead); 
}

function wordClick(e, word) { e.stopPropagation(); if(curAud){ curAud.pause(); } openAdd(word); }
function openAdd(w){ document.getElementById("mod").style.display="block"; const inp = document.getElementById("fIn"); inp.value = w; inp.removeAttribute('readonly'); fetchTranslation(w); }
function fetchTranslation(w) { fetch('/traduzir_palavra?w='+w).then(r=>r.json()).then(d=>document.getElementById('bIn').value=d.t); }
function refreshTrans() { const w = document.getElementById("fIn").value; fetchTranslation(w); }
function openManualAdd() { closeNav(); document.getElementById("mod").style.display="block"; document.getElementById("fIn").value=""; document.getElementById("fIn").removeAttribute('readonly'); document.getElementById("fIn").placeholder="Type in PT or EN (Auto-Translate)..."; }
function nextWord() { if(wordPool.length > 0) { const w = wordPool[Math.floor(Math.random() * wordPool.length)]; document.getElementById('targetWord').innerText = w.w; document.getElementById('targetIPA').innerText = w.ipa; document.getElementById('feedbackBox').innerHTML = "Tap mic & speak..."; document.getElementById('micBtn').style.background = "var(--accent)"; document.getElementById('micBtn').classList.remove('mic-btn-active'); } }

function toggleShadow(btn, e) {
    e.stopPropagation();
    if(curAud) { curAud.pause(); document.getElementById('playBtn').innerText='▶'; }
    
    if (btn.classList.contains('shadow-recording')) {
        if(shadowRec && shadowRec.state !== 'inactive') shadowRec.stop();
        btn.classList.remove('shadow-recording');
        btn.innerHTML = '🔊';
    } else {
        navigator.mediaDevices.getUserMedia({ audio: true }).then(stream => {
            shadowRec = new MediaRecorder(stream);
            shadowChunks = [];
            shadowRec.ondataavailable = e => shadowChunks.push(e.data);
            shadowRec.onstop = () => {
                const blob = new Blob(shadowChunks, { type: 'audio/ogg; codecs=opus' });
                const url = URL.createObjectURL(blob);
                const audio = new Audio(url);
                audio.play();
                btn.innerHTML = '🔊'; 
                audio.onended = () => { btn.innerHTML = '🎤'; }; 
            };
            shadowRec.start();
            btn.classList.add('shadow-recording');
            btn.innerHTML = '⏹️'; 
        }).catch(err => {
            console.error(err);
            alert("Microphone access is required for Shadowing.");
            btn.innerHTML = '🚫';
        });
    }
}

function toggleMiner() {
    document.body.classList.toggle('mode-read');
    document.body.classList.toggle('mode-mine');
    const btn = document.getElementById('modeBtn');
    if(document.body.classList.contains('mode-mine')) {
        btn.innerHTML = '⛏️ MINING MODE';
        btn.style.background = '#e67e22';
    } else {
        btn.innerHTML = '📖 READING MODE';
        btn.style.background = '#3498db';
    }
}
function mineWord(e, w) { e.stopPropagation(); openAdd(w); }
function submitAjax(e) {
    e.preventDefault();
    const f = e.target;
    const btn = f.querySelector('button');
    const originalText = btn.innerText;
    btn.innerText = "Saving...";
    fetch('/adicionar_vocab', {
        method: 'POST',
        body: new URLSearchParams(new FormData(f)) + '&ajax=1',
        headers: { 'Content-Type': 'application/x-www-form-urlencoded' }
    }).then(r=>r.json()).then(d=>{
        btn.innerText = "Saved! ✅";
        setTimeout(()=>{
            document.getElementById('mod').style.display='none';
            btn.innerText = originalText;
        }, 800);
    });
}

function loadQuiz(containerId, textContent) {
    const btn = document.getElementById('quizBtn_' + containerId);
    btn.innerText = "⏳ Generating...";
    fetch('/get_quiz', {
        method: 'POST',
        headers: {'Content-Type': 'application/x-www-form-urlencoded'},
        body: 'text=' + encodeURIComponent(textContent)
    }).then(r => r.json()).then(data => {
        const box = document.getElementById('quizBox_' + containerId);
        box.innerHTML = data.html;
        box.style.display = 'block';
        box.scrollIntoView({behavior: "smooth"});
        btn.innerText = "🔄 Regenerate Quiz";
    });
}
function reveal(el, word) {
    el.classList.remove('quiz-hidden');
    el.classList.add('quiz-revealed');
    el.innerText = word;
}
</script></head><body>
<div id="side" class="sidebar">
    <a href="javascript:void(0)" onclick="closeNav()" style="text-align:right;font-size:2rem;padding-right:20px;">&times;</a>
    <a href="javascript:void(0)" onclick="toggleTheme()">🌙 Dark Mode / Light Mode</a>
    <a href="/">📚 Play Deck <span>Estudar Cartões</span></a>
    <a href="/tutor">👩‍🏫 AI Tutor <span>Tutor Reverso</span></a>
    <a href="/vocab_list">🗂️ My Vocabulary <span>Gerenciar Palavras</span></a>
    <a href="/search">🔍 PubMed Search <span>Busca de Artigos</span></a>
    <a href="/summarizer">⚡ Smart Summarizer <span>Resumidor Inteligente</span></a>
    <a href="/pronunciation">🎤 Pronunciation <span>Treino de Fala</span></a>
    <a href="/checker">📝 Grammar & Style <span>Corretor Acadêmico</span></a>
    <a href="javascript:void(0)" onclick="openManualAdd()">➕ Add Custom Term <span>Adicionar Termo</span></a>
    <a href="/systems">🧩 Language Lab <span>Sistemas Linguísticos</span></a>
    <a href="/hub">🌐 Research Hub <span>Bases de Dados</span></a>
    <a href="/phrases">📖 Phrasebank <span>Banco de Frases</span></a>
    <a href="/miner">⛏️ Data Miner <span>Minerador de Dados</span></a>
    <a href="/novo">✨ Read & Speak <span>Leitura e Fala</span></a>
    <a href="javascript:void(0)" onclick="if(confirm('Isso apaga TUDO. Certeza?')){window.location.href='/reset_decks'}" style="color:#d63031;border-top:1px solid #eee;margin-top:20px;">⚠️ Factory Reset <span>Restaurar Padrão</span></a>
</div>
<div class="header">
    <span onclick="openNav()" style="font-size:1.5rem;cursor:pointer;">☰</span>
    <h3 style="margin-left:20px; flex:1;">{{ app_name }}</h3>
    <span id="themeIcon" onclick="toggleTheme()" style="font-size:1.5rem;cursor:pointer;">🌙</span>
</div>

{% if mode == 'list' %}
<div class="container">
    <div class="dash-row"><div class="stat-card"><div class="stat-val">{{ stats.total }}</div><div class="stat-lbl">Words</div></div><div class="stat-card"><div class="stat-val" style="color:#27ae60;">{{ stats.mastered }}</div><div class="stat-lbl">Mastered</div></div></div>
    <div class="card" style="padding:15px;"><div class="stat-lbl" style="text-align:left; margin-bottom:10px;">Activity (Last 14 days)</div><div class="heatmap">{% for day in stats.heatmap %}<div class="heat-day" style="background:{{day.color}}" title="{{day.date}}: {{day.count}}">{{day.day}}</div>{% endfor %}</div></div>
    <h3>Study Time</h3>{% for d in decks %}<a href="/jogar/{{d.id}}" style="text-decoration:none;color:inherit;"><div class="card" style="text-align:center; border:2px solid var(--accent);"><div style="font-size:3rem">{{d.icon}}</div><b>{{d.name}}</b><br><span class="sub-text">Tap to Review</span><br>{% if d.due_count > 0 %}<span class="badge" style="font-size:1rem; padding:5px 15px; margin-top:5px; display:inline-block;">{{d.due_count}} cards due</span>{% else %}<small style="color:var(--subtext)">All done!</small>{% endif %}</div></a>{% endfor %}
    <h3 style="margin-top:20px;">My Library (Saved Articles)</h3>{% for s in stories %}<div class="card" style="display:flex; justify-content:space-between; align-items:center;"><div><b>{{s.title}}</b><br><small style="color:var(--subtext)">{{s.sentences|length}} sentences</small></div><div><a href="/ler/{{s.id}}" class="btn" style="padding:5px 15px; margin-right:5px;">📖</a><a href="/deletar/{{s.id}}" class="btn" style="background:#e74c3c; padding:5px 10px;">🗑️</a></div></div>{% else %}<p style="text-align:center; color:var(--subtext);">No articles saved yet. Use 'PubMed Search' to find some!</p>{% endfor %}
</div>

{% elif mode == 'checker' %}
<div class="container"><h3>📝 Grammar & Style</h3><form action="/checker" method="POST"><div class="checker-box"><label class="checker-label">Original Text (Draft)</label><textarea name="text_input" class="text-area-box" placeholder="Paste your text here (English)..." style="height:120px;">{{original}}</textarea></div><div class="arrow-separator">⬇️ IMPROVE & POLISH ⬇️</div><button class="btn" style="width:100%; margin-bottom: 20px;">CHECK GRAMMAR</button><div class="checker-box"><label class="checker-label">Corrected Version (Academic)</label><textarea readonly id="correctedText" class="text-area-box corrected-box" style="height:120px;">{{corrected}}</textarea></div></form>{% if corrected %}<div style="margin-top:10px; display:flex; gap:10px;"><button class="btn" style="background:#34495e; flex:1;" onclick="speak(document.getElementById('correctedText').value)">🔊 LISTEN</button><button class="btn" style="background:#27ae60; flex:1;" onclick="navigator.clipboard.writeText(document.getElementById('correctedText').value); alert('Copied!')">📋 COPY</button></div>{% endif %}</div>
{% elif mode == 'tutor' %}
<div class="container"><h3>👩‍🏫 AI Tutor</h3><form action="/tutor" method="POST"><div class="checker-box"><label class="checker-label">Type in PT or EN (Auto-Detect):</label><textarea name="pt_text" class="text-area-box" placeholder="Ex: Eu preciso enviar os resultados..." style="height:100px;">{{res.pt if res else ''}}</textarea><button class="btn" style="width:100%; margin-top:10px;">TEACH ME</button></div></form>{% if res %}<div class="card" style="margin-top:20px; border-left:5px solid var(--accent);"><div class="sentence-block" style="border:none; padding:0; line-height:1.6;" onclick='prepare(this)'><span class='k-sent' style="font-size:1.2rem; font-weight:bold; cursor:pointer;">{{res.en}}</span><p style="color:var(--accent); font-weight:bold; margin-top:10px; cursor:pointer; text-align:center;">▶ TAP TO READ & LISTEN</p></div><div style="color:var(--subtext); margin-top:10px;">{{res.pt}}</div></div><h4 style="margin-top:20px;">Vocabulary Breakdown</h4><div class="card">{% for v in res.vocab %}<form action="/adicionar_vocab" method="POST" class="word-row"><div style="flex:1;"><b>{{v.word}}</b> <span style="color:var(--accent); font-family:monospace; font-size:0.8rem;">{{v.ipa}}</span><br><small style="color:var(--subtext)">{{v.trans}}</small></div><input type="hidden" name="f" value="{{v.word}}"><input type="hidden" name="b" value="{{v.trans}}"><input type="hidden" name="context" value="{{res.en}}"><div style="display:flex; gap:10px;"><button type="button" class="btn" style="padding:5px 10px; font-size:0.8rem;" onclick="speak('{{v.word}}')">🔊</button><button class="btn" style="background:#27ae60; padding:5px 10px; font-size:0.8rem;">➕</button></div></form>{% endfor %}</div>{% endif %}</div>
{% elif mode == 'vocab_list' %}
<div class="container"><h3>My Vocabulary</h3><div style="display:flex; justify-content:flex-end; gap:10px; margin-bottom:15px;"><a href="/export_vocab" class="btn" style="background:#8e44ad; text-decoration:none; font-size:0.8rem;">📥 Export</a><form action="/import_vocab" method="POST" enctype="multipart/form-data" style="display:inline;"><label for="csvUpload" class="btn" style="background:#2980b9; cursor:pointer; font-size:0.8rem;">📤 Import</label><input type="file" id="csvUpload" name="csv_file" style="display:none;" onchange="this.form.submit()"></form></div><div class="card">{% for c in cards %}<div class="word-row"><div><b>{{c.front}}</b> <span style="color:var(--accent); font-family:monospace; font-size:0.8rem;">{{c.ipa}}</span><br><small style="color:var(--subtext)">{{c.back}}</small></div><div style="display:flex; gap:10px;"><button class="btn" style="padding:5px 10px; font-size:0.8rem;" onclick="speak('{{c.front}}')">🔊</button><a href="/delete_card/{{c.id}}" class="btn" style="background:#e74c3c; padding:5px 10px; font-size:0.8rem; text-decoration:none;">🗑️</a></div></div>{% else %}<p style="text-align:center; color:var(--subtext);">No words saved yet.</p>{% endfor %}</div></div>
{% elif mode == 'search' %}
<div class="container"><h3>PubMed Search</h3><form action="/search" method="POST"><div class="checker-box"><input type="text" name="query" placeholder="e.g. 'Cancer AND Therapy'..." value="{{query}}"><input type="hidden" name="retstart" value="0"><button class="btn" style="width:100%;">SEARCH</button></div></form>{% if results %}<p style="text-align:center; color:var(--subtext); font-size:0.9rem;">Encontrados {{total_count}} resultados</p>{% for r in results %}<div class="card" style="margin-top:15px;"><b style="color:var(--accent);">{{r.journal}}</b><h4>{{r.title}}</h4><div style="max-height:80px; overflow:hidden; text-overflow:ellipsis; color:var(--subtext); font-size:0.9rem;">{{r.abstract}}</div><div style="display:flex; gap:10px; margin-top:10px;"><form action="/summarizer" method="POST" style="flex:1;"><input type="hidden" name="text_input" value="{{r.abstract}}"><input type="hidden" name="title_input" value="{{r.title}}"><button class="btn" style="background:#27ae60; width:100%; font-size:0.8rem;">⚡ PROCESS ABSTRACT</button></form><form action="/save_article" method="POST" style="flex:1;"><input type="hidden" name="abstract" value="{{r.abstract}}"><input type="hidden" name="title" value="{{r.title}}"><button class="btn" style="background:#3498db; width:100%; font-size:0.8rem;">💾 SAVE TO LIBRARY</button></form></div></div>{% endfor %}<div style="display:flex; justify-content:space-between; margin-top:20px; gap:10px;">{% if retstart > 0 %}<form action="/search" method="POST" style="flex:1;"><input type="hidden" name="query" value="{{query}}"><input type="hidden" name="retstart" value="{{retstart - 5}}"><button class="btn" style="background:#95a5a6; width:100%;">⬅ Anterior</button></form>{% endif %}{% if (retstart + 5) < total_count %}<form action="/search" method="POST" style="flex:1;"><input type="hidden" name="query" value="{{query}}"><input type="hidden" name="retstart" value="{{retstart + 5}}"><button class="btn" style="width:100%;">Próxima ➡</button></form>{% endif %}</div>{% endif %}</div>
{% elif mode == 'summarizer' %}
<div class="container"><h3>Smart Summarizer</h3><form action="/summarizer" method="POST" enctype="multipart/form-data"><div class="card"><label class="checker-label">Upload Batch (NBIB/RIS) or PDF</label><input type="file" name="nbib_file" accept=".nbib,.ris,.txt,.pdf" style="margin-bottom:10px;"><p style="color:var(--subtext); font-size:0.8rem; margin-top:0;">⚠️ Warning: Limits to first 5 articles.</p><div class="arrow-separator" style="font-size:1rem; margin:10px 0;">OR</div><label class="checker-label">Paste Text</label><textarea name="text_input" class="text-area-box" style="height:100px;" placeholder="Paste text here..."></textarea><button class="btn" style="width:100%; margin-top:15px;">PROCESS</button></div></form>{% if batch_results %}{% for res in batch_results %}<div class="card" style="margin-top:20px;"><button id="modeBtn" class="btn" style="width:100%; margin-bottom:10px; background:#3498db;" onclick="toggleMiner()">📖 READING MODE</button><h4 style="margin:0 0 10px 0;">{{res.title}}</h4><div class="sentence-block" style="border:none; padding:0; line-height:1.6;" onclick='prepare(this)'>{{ res.formatted_html | safe }}<p style="color:var(--accent); font-weight:bold; margin-top:10px; cursor:pointer; text-align:center;">▶ TAP TO READ & LISTEN</p></div><button id="quizBtn_{{loop.index}}" class="btn" style="width:100%; margin-top:15px; background:#8e44ad;" onclick='loadQuiz("{{loop.index}}", {{ res.clean_text|tojson }})'>🧩 TEST ME (Generate Quiz)</button><div id="quizBox_{{loop.index}}" style="display:none; margin-top:15px; padding:15px; background:var(--input); border:2px dashed #8e44ad; border-radius:10px; line-height:2;"></div><details style="margin-top:10px; border-top:1px solid var(--border); padding-top:10px;"><summary style="cursor:pointer; color:var(--accent); font-weight:bold;">🇧🇷 Ver Tradução</summary><p style="margin-top:10px; color:var(--text); line-height:1.5;">{{ res.translation }}</p></details></div>{% endfor %}{% endif %}</div>
{% elif mode == 'study_play' %}
<div class="container">{% if card %}<div class="card" style="text-align:center;padding:30px;"><small>{{deck.name}}</small>{% if card.is_cloze %}<h1 id="wordTxt" class="cloze-content">{{ card.cloze_hint }}</h1><p style="color:var(--subtext);font-style:italic;">Complete the sentence</p>{% else %}<h1 id="wordTxt">{{card.front}}</h1><p style="color:var(--accent);font-family:monospace;">{{card.ipa}}</p>{% endif %}<div id="ansBox" style="display:none;margin-top:20px;border-top:1px dashed #ccc;padding-top:20px;">{% if card.is_cloze %}<h2 style="color:#2ecc71;">{{card.front}}</h2><small style="display:block;margin-top:5px;color:var(--subtext);">Translation: {{card.back}}</small>{% else %}<h2 style="color:#2ecc71;">{{card.back}}</h2>{% endif %}{% if card.context %}<div style="margin-top:10px; padding:10px; background:#e1f5fe; border-radius:8px; font-style:italic; font-size:0.9rem;">💡 "{{card.context}}"</div>{% endif %}<div style="display:flex; gap:10px; margin-top:20px;"><button class="srs-btn" style="background:var(--hard);" onclick="submitRating('hard', '{{deck.id}}', '{{card.front}}')">Hard</button><button class="srs-btn" style="background:var(--med);" onclick="submitRating('medium', '{{deck.id}}', '{{card.front}}')">Good</button><button class="srs-btn" style="background:var(--easy);" onclick="submitRating('easy', '{{deck.id}}', '{{card.front}}')">Easy</button></div></div></div><div style="display:flex;flex-direction:column;gap:10px;"><button class="btn" onclick="speak('{{card.front}}')">📢 LISTEN</button><button class="btn" style="background:#34495e" onclick="document.getElementById('ansBox').style.display='block'">SHOW ANSWER</button><a href="/jogar/{{deck.id}}" class="btn" style="background:#f1c40f; color:#333; text-align:center; text-decoration:none;">⏭ Pular / Próxima</a></div>{% else %}<div class="card" style="text-align:center;padding:40px;"><h3>🎉 All caught up!</h3><p>No due cards.</p><a href="/add_random/{{deck.id}}" class="btn" style="background:#ff9f43;">🎲 ADD RANDOM</a></div>{% endif %}</div>
{% elif mode == 'pronunciation' %}
<div class="container"><h3>Pronunciation Lab</h3><div class="card" style="text-align:center; padding:40px;"><h1 id="targetWord" style="font-size:2.5rem; margin-bottom:10px; color:var(--text);">...</h1><p id="targetIPA" style="color:var(--accent); font-family:monospace; font-size:1.5rem; margin-top:0;">/ipa/</p><div id="feedbackBox" style="margin:20px 0; min-height:50px; font-weight:bold; font-size:1.1rem; color:var(--subtext);">Tap the mic and say the word.</div><button id="micBtn" class="fab" style="position:relative; bottom:auto; right:auto; width:80px; height:80px; font-size:2rem; box-shadow:0 4px 15px rgba(52, 152, 219, 0.4);" onclick="toggleRec()">🎤</button><div style="margin-top:30px; display:flex; gap:10px; justify-content:center;"><button class="btn" style="background:var(--input); color:var(--text); border:1px solid var(--border);" onclick="speak(document.getElementById('targetWord').innerText)">🔊 Listen</button><button class="btn" style="background:var(--card); color:var(--text); border:1px solid var(--border);" onclick="nextWord()">🎲 Next Word</button></div></div><script>const wordPool = {{ words_json | safe }}; const SpeechRecognition = window.SpeechRecognition || window.webkitSpeechRecognition; let recognition; function nextWord() { if(wordPool.length > 0) { const w = wordPool[Math.floor(Math.random() * wordPool.length)]; document.getElementById('targetWord').innerText = w.w; document.getElementById('targetIPA').innerText = w.ipa; document.getElementById('feedbackBox').innerHTML = "Tap mic & speak..."; document.getElementById('micBtn').style.background = "var(--accent)"; document.getElementById('micBtn').classList.remove('mic-btn-active'); } } if (SpeechRecognition) { recognition = new SpeechRecognition(); recognition.lang = 'en-US'; recognition.continuous = false; recognition.onstart = () => { document.getElementById('micBtn').classList.add('mic-btn-active'); document.getElementById('feedbackBox').innerText = "Listening..."; }; recognition.onend = () => { document.getElementById('micBtn').classList.remove('mic-btn-active'); }; recognition.onresult = (event) => { const spoken = event.results[0][0].transcript.toLowerCase().trim(); const target = document.getElementById('targetWord').innerText.toLowerCase().trim(); const confidence = Math.round(event.results[0][0].confidence * 100); const cleanSpoken = spoken.replace(/[^\\w\\s]/gi, ''); const cleanTarget = target.replace(/[^\\w\\s]/gi, ''); if (cleanSpoken === cleanTarget) { document.getElementById('feedbackBox').innerHTML = `<span style="color:#27ae60; font-size:1.4rem;">✅ Perfect!</span><br><small>Confidence: ${confidence}%</small>`; new Audio('https://actions.google.com/sounds/v1/cartoon/pop.ogg').play(); } else { document.getElementById('feedbackBox').innerHTML = `<span style="color:#e74c3c; font-size:1.4rem;">❌ Try again</span><br>You said: "<b>${spoken}</b>"`; } }; recognition.onerror = (event) => { document.getElementById('feedbackBox').innerText = "Error: " + event.error; }; } else { document.getElementById('feedbackBox').innerHTML = "⚠️ Browser not supported. Use Chrome."; } function toggleRec() { if(recognition) recognition.start(); else alert("Browser not supported"); } nextWord();</script></div>
{% elif mode == 'systems' %}
<div class="container"><h3>Language Lab</h3><form action="/systems" method="POST"><input type="text" name="word" placeholder="Type a word..."><button class="btn" style="width:100%">ANALYZE</button></form>{% if analysis %}<div class="card" style="margin-top:20px;"><div style="display:flex;justify-content:space-between;align-items:center;"><h2 style="margin:0;">{{analysis.word}}</h2><button class="btn" style="padding:10px 15px;" onclick="speak('{{analysis.word}}')">🔊</button></div><p style="color:var(--accent); font-family:monospace; font-size:1.2rem;">{{analysis.ipa}}</p><div class="sys-box"><b>🇧🇷 Translation:</b><br>{{analysis.trans}}</div><div class="sys-box"><b>🧩 Morphology:</b><br>{% for p in analysis.morph %}{{p|safe}}<br>{% else %}No common roots found.{% endfor %}</div></div>{% endif %}</div>
{% elif mode == 'hub' %}
<div class="container"><h3 style="text-align:center;">Research Hub</h3><div style="display:grid;grid-template-columns:1fr 1fr;gap:15px;">{% for link in links %}<div class="card hub-card" onclick="window.open('{{link.url}}','_blank')"><div style="font-size:2.5rem;">{{link.icon}}</div><b>{{link.name}}</b><br><small>{{link.desc}}</small></div>{% endfor %}</div></div>
{% elif mode == 'phrases' %}
<div class="container">
    <h3>Phrasebank</h3>
    <p style="font-size:0.9rem; color:var(--subtext);">Tap phrase for Audio/Karaoke. Use Mic for Shadowing.</p>
    {% for sec, list_phrases in phrases.items() %}
    <div class="card">
        <div class="section-header" style="background:#eef2f3; padding:10px; border-radius:8px; margin-bottom:10px; font-weight:bold; color:var(--accent);">{{sec}}</div>
        {% for p in list_phrases %}
        <div class="sentence-block" onclick="prepare(this)">
            <button class="btn shadow-btn" onclick="toggleShadow(this, event)">🎤</button>
            <div style="margin-right: 50px;">
                {% for w in p.en.split() %}<span class="word-span" onclick="mineWord(event, '{{w|replace("'","")|replace(".","")}}')">{{w}}</span> {% endfor %}
                <br><small style="color:var(--subtext); font-style:italic;">{{p.pt}}</small>
            </div>
        </div>
        {% endfor %}
    </div>
    {% endfor %}
</div>
{% elif mode == 'miner_home' %}
<div class="container"><form action="/miner" method="POST" enctype="multipart/form-data"><div class="card"><h3>⛏️ Data Miner</h3><label style="display:block;margin-bottom:10px;color:var(--subtext);">Suporta: PDF, TXT, RIS, NBIB</label><input type="file" name="arquivo_upload" multiple accept=".pdf,.txt,.ris,.nbib"><textarea name="texto_full" style="height:100px;margin-top:10px;" placeholder="Ou cole o texto aqui..."></textarea><button class="btn" style="width:100%">EXTRACT</button></div></form></div>
{% elif mode == 'miner_res' %}
<div class="container">{% for w, c in keywords %}<div style="display:flex;justify-content:space-between;padding:10px;border-bottom:1px solid var(--border);background:var(--card);"><span><b>{{w}}</b> ({{c}})</span><button class="btn" onclick="openAdd('{{w}}')">➕</button></div>{% endfor %}</div>
{% elif mode == 'read' %}
<div class="container">
    <div class="card" style="text-align:center; background:#f3e5f5; border:1px solid #9b59b6; margin-bottom: 20px;">
        <h3 style="margin:0 0 10px 0; color: #8e44ad;">🎧 Sci-Podcast</h3>
        <p style="font-size:0.9rem; color:#555;">Transform this text into an audio episode.</p>
        <a href="/podcast/{{story.id}}" class="btn" style="background:#8e44ad; width:100%; display:block; box-sizing:border-box; text-decoration:none;">DOWNLOAD MP3</a>
    </div>
    {% for s in story.sentences %}
    <div class="sentence-block" onclick="prepare(this)"><b>{{s.en}}</b><br><small>{{s.pt}}</small></div>
    {% endfor %}
</div>
{% elif mode == 'new' %}
<div class="container"><form action="/processar" method="POST" enctype="multipart/form-data"><div class="card"><h3>New Reading</h3><label>Upload (PDF, TXT, RIS, NBIB):</label><input type="file" name="arquivo_upload" multiple accept=".pdf,.txt,.ris,.nbib"><textarea name="texto_full" placeholder="Or paste text here..." style="height:150px; margin-top:10px;"></textarea><button class="btn" style="width:100%">PROCESS</button></div></form></div>
{% endif %}

<div class="fab" onclick="openManualAdd()">➕</div>
<div class="player"><div class="controls-row"><button id="playBtn" class="btn" style="border-radius:50%;width:50px;" onclick="togglePlay()">▶</button><select onchange="acc=this.value"><option value="com">🇺🇸</option><option value="co.uk">🇬🇧</option></select><input type="range" min="0.5" max="1.5" step="0.1" value="1.0" oninput="rate=this.value; if(curAud){curAud.playbackRate=this.value;}"></div></div>
<div id="mod" class="modal" style="display:none;"><form action="/adicionar_vocab" method="POST" onsubmit="submitAjax(event)"><h3 style="margin-top:0">Save to My Vocabulary</h3><div class="modal-actions"><input type="text" id="fIn" name="term" placeholder="Type in PT or EN (Auto-Translate)..." style="width:100%"></div><button class="btn" style="margin-top:10px;">SAVE</button><button type="button" class="btn" style="background:#888;margin-top:10px;" onclick="document.getElementById('mod').style.display='none'">CLOSE</button></form></div>
</body></html>
"""

if __name__ == '__main__': app.run(host='0.0.0.0', port=5000)
