import re
import json
import asyncio
import edge_tts
from io import BytesIO
from collections import Counter
from deep_translator import GoogleTranslator
import eng_to_ipa as ipa

# --- CONFIGURAÇÕES DE VOZ (NEURAL) ---
VOICE_MAPPING = {
    'com': 'en-US-ChristopherNeural',  # Voz Americana Acadêmica
    'co.uk': 'en-GB-RyanNeural',       # Voz Britânica
    'pt': 'pt-BR-AntonioNeural'        # Voz Português (Opcional)
}

# --- BANCO DE DADOS DE PESQUISA ---
RESEARCH_LINKS = [
    {"name": "PubMed", "url": "https://pubmed.ncbi.nlm.nih.gov/", "icon": "🧬"},
    {"name": "SciELO", "url": "https://scielo.org/", "icon": "🌎"},
    {"name": "Google Scholar", "url": "https://scholar.google.com.br/", "icon": "🎓"},
    {"name": "Connected Papers", "url": "https://www.connectedpapers.com/", "icon": "🕸️"},
    {"name": "Semantic Scholar", "url": "https://www.semanticscholar.org/", "icon": "🧠"}
]

# --- PHRASEBANK PARA O SCI-WRITER ---
ACADEMIC_PHRASEBANK = {
    "1. Introduction": [
        {"en": "The primary objective of this study is to...", "pt": "O objetivo principal é..."},
        {"en": "Recent developments have heightened the need for...", "pt": "Desenvolvimentos recentes aumentaram a necessidade de..."},
        {"en": "This paper addresses the issue of...", "pt": "Este artigo aborda a questão de..."}
    ],
    "2. Contrast/Argue": [
        {"en": "However, this approach has limitations.", "pt": "No entanto, tem limitações."},
        {"en": "On the other hand,", "pt": "Por outro lado,"},
        {"en": "Conversely,", "pt": "Inversamente,"},
        {"en": "While preliminary results suggest..., it remains unclear...", "pt": "Embora sugiram..., não está claro..."}
    ],
    "3. Addition": [
        {"en": "Furthermore,", "pt": "Além disso,"},
        {"en": "In addition to...", "pt": "Em adição a..."},
        {"en": "Moreover,", "pt": "Além do mais,"}
    ],
    "4. Methodology": [
        {"en": "Data were collected using...", "pt": "Dados coletados usando..."},
        {"en": "The participants were divided into...", "pt": "Participantes divididos em..."},
        {"en": "Statistical analysis was performed using...", "pt": "Análise feita com..."}
    ],
    "5. Results/Conclusion": [
        {"en": "The results indicate that...", "pt": "Os resultados indicam que..."},
        {"en": "These findings highlight the importance of...", "pt": "Estes achados destacam a importância de..."},
        {"en": "In conclusion,", "pt": "Em conclusão,"}
    ]
}

# --- FUNÇÕES DE ÁUDIO ---
async def generate_neural_audio(text, voice):
    """Gera áudio usando a API Edge TTS (Microsoft)"""
    communicate = edge_tts.Communicate(text, voice)
    audio_data = BytesIO()
    async for chunk in communicate.stream():
        if chunk["type"] == "audio":
            audio_data.write(chunk["data"])
    audio_data.seek(0)
    return audio_data

def get_audio_sync(text, accent='com'):
    """Wrapper síncrono para chamar o async no Flask"""
    if not text: return None
    voice = VOICE_MAPPING.get(accent, 'en-US-ChristopherNeural')
    try:
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        return loop.run_until_complete(generate_neural_audio(text, voice))
    except Exception as e:
        print(f"Erro TTS: {e}")
        return None

# --- FERRAMENTAS DE TEXTO ---
def format_abstract_smart(text):
    """
    Formata texto cru de PDF para HTML legível.
    1. Remove quebras de linha estranhas.
    2. Coloca em negrito seções comuns (Introduction, Methods...).
    3. Separa em blocos de parágrafos.
    """
    if not text: return ""
    
    # Limpa espaços múltiplos
    formatted = re.sub(r'\s+', ' ', text)
    
    # Lista de seções para destacar
    sections = [
        "BACKGROUND", "OBJECTIVE", "METHODS", "RESULTS", "CONCLUSIONS", 
        "CONCLUSION", "DISCUSSION", "INTRODUCTION", "REFERENCES"
    ]
    
    # Adiciona quebra de linha e negrito antes das seções
    for sec in sections:
        # Regex procura a palavra seguida de : ou espaço, case insensitive
        pattern = re.compile(rf"({sec}[:\s])", re.IGNORECASE)
        formatted = pattern.sub(r"<br><br><b style='color:#2c3e50; font-size:1.1em;'>\1</b>", formatted)
    
    return formatted

def get_phonetic(text):
    """Gera IPA para uma palavra ou frase"""
    if not text: return ""
    # Remove pontuação para o conversor não quebrar
    clean = re.sub(r'[^\w\s]', '', text)
    return "/" + ipa.convert(clean) + "/"

# --- HELPERS DIVERSOS ---
def get_top_keywords(text):
    """Extrai palavras-chave para o Minerador"""
    words = re.findall(r'\b[a-zA-Z]{5,}\b', text.lower())
    stop_words = {"which", "their", "about", "would", "these", "other", "words", "could", "write", "first", "water", "after", "where", "right", "think", "three", "years", "place", "sound", "great", "again", "still", "study", "using", "group", "results", "table", "level", "based", "found", "value", "total", "during", "between", "analysis"}
    filtered = [w for w in words if w not in stop_words]
    return Counter(filtered).most_common(20)
