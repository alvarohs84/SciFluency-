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

# Banco de Frases
ACADEMIC_PHRASEBANK = {
    "1. Introduction": [{"en": "Recent developments in this field have heightened the need for...", "pt": "Desenvolvimentos recentes..."}],
    "2. Methods": [{"en": "Data were collected using a semi-structured interview...", "pt": "Dados coletados via..."}],
    "3. Results": [{"en": "The results indicate that...", "pt": "Os resultados indicam que..."}],
    "4. Discussion": [{"en": "These findings suggest that...", "pt": "Esses achados sugerem..."}]
}

ACADEMIC_REPLACEMENTS = {
    "big": "substantial", "huge": "significant", "bad": "detrimental",
    "good": "beneficial", "think": "hypothesize", "get": "obtain",
    "make": "generate", "show": "demonstrate"
}

# --- PARSER DE ARQUIVOS BIBLIOGRÁFICOS (RIS / NBIB) ---
def parse_bib_file(content):
    """Lê texto RIS/NBIB e retorna lista de referências"""
    references = []
    current = {}
    
    # Normaliza quebras de linha
    lines = content.replace('\r\n', '\n').split('\n')
    
    for line in lines:
        line = line.strip()
        if not line: continue
        
        # Detecta fim de registro (ER no RIS)
        if line.startswith("ER  -"):
            if current.get('title'): references.append(current)
            current = {}
            continue
            
        # Título (TI, T1)
        if line.startswith("TI  - ") or line.startswith("T1  - "):
            current['title'] = line[6:]
            
        # Autor (AU, A1, FAU) - Pega o primeiro
        elif (line.startswith("AU  - ") or line.startswith("A1  - ") or line.startswith("FAU - ")) and 'authors' not in current:
            current['authors'] = line[6:]
            
        # Ano (PY, Y1, DP)
        elif (line.startswith("PY  - ") or line.startswith("Y1  - ") or line.startswith("DP  - ")):
            # Tenta extrair 4 digitos
            match = re.search(r'\d{4}', line)
            if match: current['year'] = match.group(0)
            
        # Resumo (AB, N2)
        elif line.startswith("AB  - ") or line.startswith("N2  - "):
            current['abstract'] = line[6:]

    # Adiciona o último se não tiver fechado
    if current.get('title'): references.append(current)
    return references

# --- CITAÇÕES ---
def generate_citation_formats(ref):
    if not ref: return {}
    aut = ref.authors if ref.authors else "AUTOR"
    title = ref.title if ref.title else "Título"
    year = ref.year if ref.year else "s.d."
    return {
        "abnt": f"{aut.upper()}. <b>{title}</b>. {year}.",
        "apa": f"{aut}. ({year}). <i>{title}</i>.",
        "vancouver": f"{aut}. {title}. {year}.",
        "ama": f"{aut}. {title}. {year}.",
        "ieee": f"[1] {aut}, \"{title},\" {year}."
    }

# --- AUXILIARES (TTS, FORMAT) ---
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
    formatted = re.sub(r'\s+', ' ', text).strip()
    for sec in ["BACKGROUND", "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSIONS", "INTRODUCTION"]:
        formatted = re.sub(rf"({sec}[:\s])", r"<br><br><b style='color:#2980b9'>\1</b>", formatted, flags=re.IGNORECASE)
    return formatted

def improve_english_text(text):
    fixed = text
    for s, a in ACADEMIC_REPLACEMENTS.items():
        fixed = re.sub(rf"\b{s}\b", f"<b>{a}</b>", fixed, flags=re.IGNORECASE)
    return fixed

def get_top_keywords(text):
    words = re.findall(r'\b[a-zA-Z]{5,}\b', text.lower())
    stop = {"which", "about", "their", "these", "other", "after", "where", "study", "using", "results", "group"}
    return Counter([w for w in words if w not in stop]).most_common(20)

def search_pubmed(query):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
        ids = Entrez.read(handle)['IdList']
        if not ids: return []
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        raw = handle.read()
        return parse_bib_file(raw) # Reusa o parser criado acima!
    except: return []

def get_phonetic(text):
    if not text: return ""
    clean = re.sub(r'[^\w\s]', '', text)
    try: return "/" + ipa.convert(clean) + "/"
    except: return ""
