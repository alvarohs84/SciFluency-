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

# --- GERADOR DE CITAÇÕES (8 ESTILOS) ---
def generate_citation_formats(ref):
    """Gera strings de citação (ABNT, APA, Vancouver, AMA, NLM, Harvard, IEEE, Chicago)"""
    if not ref: return {}
    
    # Dados brutos
    aut = ref.authors if ref.authors and ref.authors != "Unknown" else "AUTOR, A."
    title = ref.title if ref.title else "Título do Artigo"
    year = ref.year if ref.year else "s.d."
    
    # 1. ABNT (Brasil - Nomes em maiúsculo)
    abnt = f"{aut.upper()}. <b>{title}</b>. {year}."
    
    # 2. APA (EUA - Título em itálico)
    apa = f"{aut}. ({year}). <i>{title}</i>."
    
    # 3. Vancouver (Biomédica Geral)
    vancouver = f"{aut}. {title}. {year}."

    # 4. AMA (American Medical Association)
    ama = f"{aut}. {title}. {year}."

    # 5. NLM (National Library of Medicine)
    nlm = f"{aut}. {title}. {year}."

    # 6. Harvard (Autor-Data)
    harvard = f"{aut} ({year}) '{title}'."

    # 7. IEEE (Engenharia)
    ieee = f"[1] {aut}, \"{title},\" {year}."

    # 8. Chicago (Humanidades)
    chicago = f"{aut}. \"{title}.\" {year}."
    
    return {
        "abnt": abnt, "apa": apa, "vancouver": vancouver,
        "ama": ama, "nlm": nlm, "harvard": harvard,
        "ieee": ieee, "chicago": chicago
    }

# --- FUNÇÕES DE ÁUDIO E TEXTO ---

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
    for sec in ["BACKGROUND", "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSIONS", "DISCUSSION", "INTRODUCTION"]:
        formatted = re.sub(rf"({sec}[:\s])", r"<br><br><b style='color:#2980b9'>\1</b>", formatted, flags=re.IGNORECASE)
    return formatted

def improve_english_text(text):
    fixed = text
    for s, a in ACADEMIC_REPLACEMENTS.items():
        fixed = re.sub(rf"\b{s}\b", f"<b>{a}</b>", fixed, flags=re.IGNORECASE)
    return fixed

def get_top_keywords(text):
    words = re.findall(r'\b[a-zA-Z]{5,}\b', text.lower())
    stop = {"which", "about", "their", "these", "other", "after", "where", "would", "could", "study", "using", "results", "group", "based", "table", "value"}
    return Counter([w for w in words if w not in stop]).most_common(20)

def search_pubmed(query):
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
        record = Entrez.read(handle)
        ids = record['IdList']
        if not ids: return []
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        raw = handle.read()
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

def get_phonetic(text):
    if not text: return ""
    clean = re.sub(r'[^\w\s]', '', text)
    try: return "/" + ipa.convert(clean) + "/"
    except: return ""
