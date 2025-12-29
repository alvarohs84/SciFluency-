# ... (código anterior mantido)

def get_phonetic(text):
    """Gera IPA para uma palavra ou frase"""
    if not text: return ""
    # Remove pontuação para o conversor não quebrar
    clean = re.sub(r'[^\w\s]', '', text)
    try:
        # Tenta converter. Se a lib retornar *, usa a palavra original
        return "/" + ipa.convert(clean) + "/"
    except:
        return ""
