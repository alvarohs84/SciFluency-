import re
import random
from deep_translator import GoogleTranslator
import eng_to_ipa as ipa
from collections import Counter
from pypdf import PdfReader

# Constantes e Dicionários
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
            clean_w = re.sub(r'[^\w]', '', w)
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
                clean_w = re.sub(r"[^\w']", "", word)
                if clean_w:
                    words_html += f"<span class='word-span' onclick='mineWord(event, \"{clean_w}\")'>{word}</span> "
                else:
                    words_html += word + " "
            final_html += f"<span class='k-sent' onclick='prepare(this, \"{sent}\")'>{words_html}</span>"
    return final_html
