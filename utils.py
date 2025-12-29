import re
import json
import asyncio
import edge_tts
from io import BytesIO
from collections import Counter
from deep_translator import GoogleTranslator
from pypdf import PdfReader
from Bio import Entrez
import eng_to_ipa as ipa

# --- BANCO DE FRASES ACADÊMICAS (EXPANDIDO) ---
ACADEMIC_PHRASEBANK = {
    "1. Start": [
        {"en": "The primary objective of this study is to...", "pt": "O objetivo principal deste estudo é..."},
        {"en": "Recent developments have heightened the need for...", "pt": "Desenvolvimentos recentes aumentaram a necessidade de..."},
        {"en": "Little is known about...", "pt": "Pouco se sabe sobre..."},
        {"en": "This paper addresses the issue of...", "pt": "Este artigo aborda a questão de..."}
    ],
    "2. Contrast": [
        {"en": "However, this approach has limitations.", "pt": "No entanto, esta abordagem tem limitações."},
        {"en": "On the other hand,", "pt": "Por outro lado,"},
        {"en": "Conversely,", "pt": "Inversamente,"},
        {"en": "In contrast to earlier findings,", "pt": "Em contraste com achados anteriores,"},
        {"en": "While preliminary results suggest..., it is unclear...", "pt": "Embora resultados preliminares sugiram..., não está claro..."}
    ],
    "3. Add": [
        {"en": "Furthermore,", "pt": "Além disso,"},
        {"en": "Moreover,", "pt": "Além do mais,"},
        {"en": "In addition to...", "pt": "Em adição a..."},
        {"en": "Similarly,", "pt": "Similarmente,"},
        {"en": "Another key factor is...", "pt": "Outro fator chave é..."}
    ],
    "4. Cause": [
        {"en": "Therefore,", "pt": "Portanto,"},
        {"en": "Consequently,", "pt": "Consequentemente,"},
        {"en": "As a result,", "pt": "Como resultado,"},
        {"en": "This suggests that...", "pt": "Isso sugere que..."}
    ],
    "5. Methods": [
        {"en": "Data were collected using...", "pt": "Os dados foram coletados usando..."},
        {"en": "The participants were divided into...", "pt": "Os participantes foram divididos em..."},
        {"en": "Statistical analysis was performed using...", "pt": "A análise estatística foi feita usando..."}
    ],
    "6. End": [
        {"en": "In conclusion,", "pt": "Em conclusão,"},
        {"en": "These findings highlight the importance of...", "pt": "Esses achados destacam a importância de..."},
        {"en": "Future research should focus on...", "pt": "Pesquisas futuras devem focar em..."}
    ]
}

ACADEMIC_REPLACEMENTS = {
    "big": "substantial", "huge": "significant", "bad": "detrimental",
    "good": "beneficial", "think": "hypothesize", "get": "obtain", "use": "utilize",
    "look at": "examine", "show": "demonstrate", "prove": "validate"
}
MORPHOLOGY_DB = {"un": "Not (Não)", "re": "Again (Novamente)", "itis": "Inflammation", "logy": "Study of"}

# --- FUNÇÕES TÉCNICAS ---
VOICE_MAPPING = {'com': 'en-US-ChristopherNeural', 'co.uk': 'en-GB-RyanNeural', 'pt': 'pt-BR-AntonioNeural'}

async def generate_neural_audio(text, voice):
    communicate = edge_tts.Communicate(text, voice)
    audio_data = BytesIO()
    async for chunk in communicate.stream():
        if chunk["type"] == "audio": audio_data.write(chunk["data"])
    audio_data.seek(0)
    return audio_data

def get_audio_sync(text, accent='com'):
    voice = VOICE_MAPPING.get(accent, 'en-US-ChristopherNeural')
    try:
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        return loop.run_until_complete(generate_neural_audio(text, voice))
    except: return None

def get_phonetic(text):
    if not text: return ""
    return "/" + " ".join([ipa.convert(re.sub(r'[^\w\s]', '', w)) if "*" not in ipa.convert(re.sub(r'[^\w\s]', '', w)) else w for w in text.split()]) + "/"

def process_language_logic(text):
    clean = text.strip()
    try:
        en = GoogleTranslator(source='pt', target='en').translate(clean)
        if en.lower() != clean.lower(): return en, clean
        pt = GoogleTranslator(source='en', target='pt').translate(clean)
        return clean, pt
    except: return clean, "Error"

def improve_english_text(text):
    try:
        fixed = text
        for s, a in ACADEMIC_REPLACEMENTS.items():
            fixed = re.sub(rf"\b{s}\b", a, fixed, flags=re.IGNORECASE)
        return fixed
    except: return text

def generate_smart_quiz(text):
    words = text.split()
    output = []
    for i, w in enumerate(words):
        if len(w) > 4 and i % 4 == 0: output.append(f"<span class='quiz-hidden' onclick='reveal(this, \"{w}\")'>[ ? ]</span>")
        else: output.append(w)
    return " ".join(output)

def parse_ris_nbib(content):
    return re.sub(r'\n', ' ', content)

def get_top_keywords(text):
    words = re.findall(r'\b[a-zA-Z]{5,}\b', text.lower())
    return Counter(words).most_common(20)

def parse_nbib_bulk(content):
    # Lógica de parser simplificada
    articles = []
    lines = content.split('\n')
    curr = {"title":"", "abstract":""}
    for l in lines:
        if l.startswith("TI  - "): curr["title"] = l[6:]
        if l.startswith("AB  - "): curr["abstract"] += l[6:]
        if l.startswith("PMID-"): 
            if curr["title"]: articles.append(curr)
            curr = {"title":"", "abstract":""}
    return articles[:10]

def format_abstract_smart(text):
    formatted = re.sub(r'\s+', ' ', text)
    for sec in ["BACKGROUND", "METHODS", "RESULTS", "CONCLUSIONS"]:
        formatted = re.sub(rf"({sec}[:\s])", r"<br><b>\1</b>", formatted, flags=re.IGNORECASE)
    sentences = re.split(r'(?<=[.!?])\s+', formatted)
    # Gera blocos de frase limpos, sem spans individuais para evitar poluição em PDF
    return "".join([f"<div class='sentence-block' onclick='prepare(this)' style='margin-bottom:8px;padding:5px;cursor:pointer;'>{s}</div>" for s in sentences if s])
