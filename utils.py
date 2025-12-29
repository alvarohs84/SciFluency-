import re
import json
import asyncio
import edge_tts
from io import BytesIO
from deep_translator import GoogleTranslator
import eng_to_ipa as ipa

# Vozes Humanas (Mantidas pois são melhores)
VOICE_MAPPING = {'com': 'en-US-ChristopherNeural', 'co.uk': 'en-GB-RyanNeural'}

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

def format_text_smart(text):
    # Formata texto para leitura agradável
    text = re.sub(r'\s+', ' ', text).strip()
    sentences = re.split(r'(?<=[.!?])\s+', text)
    html = ""
    for s in sentences:
        if len(s) > 2:
            # Cria span para cada palavra para permitir o clique individual
            words_html = " ".join([f"<span class='word-span' onclick='mineWord(event, \"{re.sub(r'[^\w]', '', w)}\")'>{w}</span>" for w in s.split()])
            html += f"<div class='sentence-block' onclick='prepare(this)'>{words_html}</div>"
    return html

def get_phonetic(text):
    clean = re.sub(r'[^\w\s]', '', text)
    return "/" + ipa.convert(clean) + "/"
