import re
import json
import asyncio
import edge_tts
from io import BytesIO
from collections import Counter
from deep_translator import GoogleTranslator
from pypdf import PdfReader
from gtts import gTTS
from Bio import Entrez
import eng_to_ipa as ipa
import time

# --- CONFIGURAÇÕES ---
Entrez.email = "researcher@example.com"
Entrez.tool = "SciFluencyResearch"

# --- DADOS CONSTANTES ---
VOICE_MAPPING = {
    'com': 'en-US-ChristopherNeural',
    'co.uk': 'en-GB-RyanNeural',
    'pt': 'pt-BR-AntonioNeural'
}

ACADEMIC_PHRASEBANK = {
    "1. Introduction & Context": [
        {"en": "Recent developments in this field have heightened the need for...", "pt": "Desenvolvimentos recentes neste campo aumentaram a necessidade de..."},
        {"en": "Currently, there is a paucity of data regarding...", "pt": "Atualmente, há escassez de dados sobre..."},
        {"en": "This study aims to investigate the relationship between...", "pt": "Este estudo visa investigar a relação entre..."},
        {"en": "Previous research has established that...", "pt": "Pesquisas anteriores estabeleceram que..."}
    ],
    "2. Methods & Materials": [
        {"en": "Data were collected using a semi-structured interview guide.", "pt": "Os dados foram coletados usando um roteiro de entrevista semiestruturado."},
        {"en": "The participants were divided into two groups.", "pt": "Os participantes foram divididos em dois grupos."},
        {"en": "Statistical analysis was performed using SPSS software.", "pt": "A análise estatística foi realizada usando o software SPSS."}
    ],
    "3. Results & Findings": [
        {"en": "There was a significant correlation between...", "pt": "Houve uma correlação significativa entre..."},
        {"en": "The results indicate that...", "pt": "Os resultados indicam que..."}
    ],
    "4. Discussion & Argumentation": [
        {"en": "These findings suggest that...", "pt": "Esses achados sugerem que..."},
        {"en": "However, some limitations should be noted.", "pt": "No entanto, algumas limitações devem ser notadas."}
    ],
    "5. Conclusion": [
        {"en": "In conclusion, this study demonstrates that...", "pt": "Em conclusão, este estudo demonstra que..."},
        {"en": "The evidence from this study suggests...", "pt": "As evidências deste estudo sugerem..."}
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

# --- FUNÇÕES ---
async def generate_neural_audio(text, voice):
    communicate = edge_tts.Communicate(text, voice)
    audio_data = BytesIO()
    async for chunk in communicate.stream():
        if chunk["type"] == "audio":
            audio_data.write(chunk["data"])
    audio_data.seek(0)
    return audio_data

def get_audio_sync(text, accent='com'):
    voice = VOICE_MAPPING.get(accent, 'en-US-ChristopherNeural')
    try:
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        return loop.run_until_complete(generate_neural_audio(text, voice))
    except Exception as e:
        print(f"Erro Neural TTS: {e}")
        return None

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
    formatted = re.sub(r'\s+', ' ', formatted)
    sections = ["BACKGROUND", "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSIONS", "CONCLUSION", "DISCUSSION"]
    for sec in sections:
        pattern = re.compile(rf"({sec}[:\s])", re.IGNORECASE)
        formatted = pattern.sub(r"<br><br><b>\1</b>", formatted)
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
            final_html += f"<div class='sentence-block' onclick='prepare(this)' style='margin-bottom:8px; padding:8px; cursor:pointer; line-height:1.6;'>{words_html}</div>"
    return final_html
