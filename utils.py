import re
import json
import asyncio
import edge_tts
from io import BytesIO
from collections import Counter
from deep_translator import GoogleTranslator
from Bio import Entrez
import eng_to_ipa as ipa

# Configurações
Entrez.email = "student@scifluency.com"
VOICE_MAPPING = {'com': 'en-US-ChristopherNeural', 'co.uk': 'en-GB-RyanNeural'}

# Banco de Frases (Writer)
ACADEMIC_PHRASEBANK = {
    "1. Introduction": [{"en": "Recent developments in this field have heightened the need for...", "pt": "Desenvolvimentos recentes..."}],
    "2. Methods": [{"en": "Data were collected using a semi-structured interview...", "pt": "Dados coletados via..."}],
    "3. Results": [{"en": "The results indicate that...", "pt": "Os resultados indicam que..."}],
    "4. Discussion": [{"en": "These findings suggest that...", "pt": "Esses achados sugerem..."}],
    "5. Connectors": [{"en": "Furthermore,", "pt": "Além disso,"}, {"en": "However,", "pt": "No entanto,"}, {"en": "Therefore,", "pt": "Portanto,"}]
}

# Banco de Substituições (Grammar Checker)
ACADEMIC_REPLACEMENTS = {
    "big": "substantial", "huge": "significant", "bad": "detrimental",
    "good": "beneficial", "think": "hypothesize", "get": "obtain",
    "make": "generate", "show": "demonstrate", "use": "utilize",
    "really": "significantly", "very": "highly", "look at": "examine"
}

# --- FUNÇÕES ---

async def generate_neural_audio(text, voice):
    communicate = edge_tts.Communicate(text, voice)
    audio_data = BytesIO()
    async for chunk in communicate.stream():
        if chunk["type"] == "audio": audio_data.write(chunk["data"])
    audio_data.seek(0)
    return audio_data

def get_audio_sync(text, accent='com'):
    if not text: return None
    voice = VOICE_MAPPING.get(accent, 'en-US-ChristopherNeural')
    try:
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        return loop.run_until_complete(generate_neural_audio(text, voice))
    except: return None

def format_abstract_smart(text):
    # Formata PDF cru para HTML bonito com negrito nas seções
    formatted = re.sub(r'\s+', ' ', text).strip()
    for sec in ["BACKGROUND", "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSIONS", "DISCUSSION", "INTRODUCTION"]:
        formatted = re.sub(rf"({sec}[:\s])", r"<br><br><b style='color:#2980b9'>\1</b>", formatted, flags=re.IGNORECASE)
    
    # Quebra em blocos de frase para o Karaokê
    sentences = re.split(r'(?<=[.!?])\s+', formatted)
    html = ""
    for s in sentences:
        if len(s) > 2:
            # Spans para cada palavra (Mineração)
            words_html = " ".join([f"<span class='word-span' onclick='mineWord(event, \"{re.sub(r'[^\w]', '', w)}\")'>{w}</span>" for w in s.split()])
            html += f"<div class='sentence-block' onclick='prepare(this)'>{words_html}</div>"
    return html

def improve_english_text(text):
    # Grammar Checker Básico
    fixed = text
    for s, a in ACADEMIC_REPLACEMENTS.items():
        fixed = re.sub(rf"\b{s}\b", f"<b>{a}</b>", fixed, flags=re.IGNORECASE)
    return fixed

def get_top_keywords(text):
    # Data Miner
    words = re.findall(r'\b[a-zA-Z]{5,}\b', text.lower())
    stop = {"which", "about", "their", "these", "other", "after", "where", "would", "could", "study", "using", "results", "group", "based", "table", "value"}
    return Counter([w for w in words if w not in stop]).most_common(20)

def search_pubmed(query):
    # Busca real no PubMed
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
        record = Entrez.read(handle)
        ids = record['IdList']
        if not ids: return []
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        raw = handle.read()
        # Parser simples de Medline
        articles = []
        curr = {"title":"", "abstract":"", "journal": "PubMed"}
        for line in raw.split('\n'):
            if line.startswith("TI  - "): curr['title'] = line[6:]
            elif line.startswith("AB  - "): curr['abstract'] += line[6:]
            elif line.startswith("PMID-"): 
                if curr['title']: articles.append(curr)
                curr = {"title":"", "abstract":"", "journal": "PubMed"}
        if curr['title']: articles.append(curr)
        return articles
    except: return []
