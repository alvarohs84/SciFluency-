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

# Frases e Substituições
ACADEMIC_PHRASEBANK = {
    "1. Introduction": [{"en": "Recent developments...", "pt": "Desenvolvimentos..."}],
    "2. Methods": [{"en": "Data were collected...", "pt": "Dados coletados..."}],
    "3. Results": [{"en": "The results indicate...", "pt": "Resultados indicam..."}],
    "4. Discussion": [{"en": "These findings suggest...", "pt": "Achados sugerem..."}]
}
ACADEMIC_REPLACEMENTS = {"big": "substantial", "huge": "significant", "bad": "detrimental", "good": "beneficial", "think": "hypothesize"}

# --- PARSER RIS/NBIB ---
def parse_bib_file(content):
    references = []
    current = {}
    lines = content.replace('\r\n', '\n').split('\n')
    for line in lines:
        line = line.strip()
        if not line: continue
        if line.startswith("ER  -"):
            if current.get('title'): references.append(current)
            current = {}
            continue
        if line.startswith("TI  - ") or line.startswith("T1  - "): current['title'] = line[6:]
        elif (line.startswith("AU  - ") or line.startswith("A1  - ")) and 'authors' not in current: current['authors'] = line[6:]
        elif (line.startswith("PY  - ") or line.startswith("DP  - ")):
            match = re.search(r'\d{4}', line)
            if match: current['year'] = match.group(0)
        elif line.startswith("AB  - ") or line.startswith("N2  - "): current['abstract'] = line[6:]
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
        "nlm": f"{aut}. {title}. {year}.",
        "harvard": f"{aut} ({year}) '{title}'.",
        "ieee": f"[1] {aut}, \"{title},\" {year}.",
        "chicago": f"{aut}. \"{title}.\" {year}."
    }

# --- PUBMED COM PAGINAÇÃO ---
def search_pubmed(query, start=0):
    if not query: return []
    try:
        # retstart define o início da lista
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10, retstart=start, sort="relevance")
        id_list = Entrez.read(handle)['IdList']
        handle.close()
        if not id_list: return []
        
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
        raw_data = handle.read()
        handle.close()
        
        articles = []
        curr = {}
        for line in raw_data.strip().split('\n'):
            key = line[:4].strip()
            val = line[6:].strip()
            if key == 'TI': curr['title'] = val
            elif key == '' and 'title' in curr and 'abstract' not in curr: curr['title'] += " " + val
            elif key == 'AB': curr['abstract'] = val
            elif key == '' and 'abstract' in curr: curr['abstract'] += " " + val
            elif key == 'TA': curr['journal'] = val
            elif key == 'PMID':
                if 'title' in curr:
                    curr['url'] = f"https://pubmed.ncbi.nlm.nih.gov/{curr['pmid']}/"
                    articles.append(curr)
                curr = {'pmid': val, 'title': '...', 'abstract': '...', 'journal': 'PubMed'}
        if 'title' in curr:
             curr['url'] = f"https://pubmed.ncbi.nlm.nih.gov/{curr['pmid']}/"
             articles.append(curr)
        return articles
    except Exception as e:
        print(e)
        return []

# --- AUXILIARES ---
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
    for sec in ["BACKGROUND", "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSIONS", "DISCUSSION"]:
        formatted = re.sub(rf"({sec}[:\s])", r"<br><br><b style='color:#2980b9'>\1</b>", formatted, flags=re.IGNORECASE)
    return formatted

def improve_english_text(text):
    fixed = text
    for s, a in ACADEMIC_REPLACEMENTS.items(): fixed = re.sub(rf"\b{s}\b", f"<b>{a}</b>", fixed, flags=re.IGNORECASE)
    return fixed

def get_top_keywords(text):
    words = re.findall(r'\b[a-zA-Z]{5,}\b', text.lower())
    stop = {"which", "about", "their", "these", "other", "after", "where", "study", "using", "results", "group"}
    return Counter([w for w in words if w not in stop]).most_common(20)

def get_phonetic(text):
    if not text: return ""
    clean = re.sub(r'[^\w\s]', '', text)
    try: return "/" + ipa.convert(clean) + "/"
    except: return ""
