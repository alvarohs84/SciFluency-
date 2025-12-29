from flask_sqlalchemy import SQLAlchemy
from datetime import datetime

db = SQLAlchemy()

# --- SISTEMA DE APRENDIZADO (SRS) ---
class Deck(db.Model):
    id = db.Column(db.String(50), primary_key=True)
    name = db.Column(db.String(100))
    pt_name = db.Column(db.String(100))
    icon = db.Column(db.String(10))
    cards = db.relationship('Card', backref='deck', lazy=True, cascade="all, delete-orphan")

class Card(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    front = db.Column(db.String(200)) 
    back = db.Column(db.String(200))
    ipa = db.Column(db.String(200))
    context = db.Column(db.Text)
    deck_id = db.Column(db.String(50), db.ForeignKey('deck.id'))
    
    # SRS Math
    next_review = db.Column(db.String(20)) # YYYY-MM-DD
    interval = db.Column(db.Integer, default=0)
    ease_factor = db.Column(db.Float, default=2.5) # Para algoritmo SM-2 futuro

# --- SISTEMA DE PESQUISA (THESIS OS) ---
class Project(db.Model):
    """Um TCC, Dissertação ou Artigo sendo escrito"""
    id = db.Column(db.String(50), primary_key=True) # ex: 'mestrado'
    title = db.Column(db.String(200))
    target_journal = db.Column(db.String(100)) # Onde pretende publicar
    references = db.relationship('Reference', backref='project', lazy=True)
    drafts = db.relationship('Draft', backref='project', lazy=True)

class Reference(db.Model):
    """Um artigo/livro salvo na biblioteca"""
    id = db.Column(db.Integer, primary_key=True)
    title = db.Column(db.String(300))
    authors = db.Column(db.String(200))
    year = db.Column(db.String(4))
    journal = db.Column(db.String(100))
    abstract = db.Column(db.Text)
    
    # Arquivo e Status
    pdf_filename = db.Column(db.String(200)) # Nome do arquivo salvo
    status = db.Column(db.String(20)) # 'to_read', 'reading', 'done'
    
    # Relacionamentos
    project_id = db.Column(db.String(50), db.ForeignKey('project.id'))
    notes = db.relationship('Note', backref='reference', lazy=True, cascade="all, delete-orphan")

class Note(db.Model):
    """Fichamento: Uma anotação sobre um trecho específico"""
    id = db.Column(db.Integer, primary_key=True)
    content = db.Column(db.Text) # O que o aluno escreveu
    quote = db.Column(db.Text)   # O trecho original do PDF (citação direta)
    page_num = db.Column(db.Integer)
    tags = db.Column(db.String(100)) # ex: 'metodologia', 'importante'
    reference_id = db.Column(db.Integer, db.ForeignKey('reference.id'))

class Draft(db.Model):
    """Os textos escritos pelo aluno (Intro, Methods, etc)"""
    id = db.Column(db.Integer, primary_key=True)
    section_name = db.Column(db.String(50)) # Introduction, Methods...
    content = db.Column(db.Text) # O texto em si (HTML/Markdown)
    last_updated = db.Column(db.DateTime, default=datetime.utcnow)
    project_id = db.Column(db.String(50), db.ForeignKey('project.id'))

# --- LEGADO E UTILITÁRIOS ---
class Story(db.Model):
    # Mantido para compatibilidade com o "Podcast Mode" atual
    id = db.Column(db.String(50), primary_key=True)
    title = db.Column(db.String(200))
    sentences = db.relationship('Sentence', backref='story', lazy=True, cascade="all, delete-orphan")

class Sentence(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    en = db.Column(db.Text)
    pt = db.Column(db.Text)
    story_id = db.Column(db.String(50), db.ForeignKey('story.id'))

class StudyLog(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    date = db.Column(db.String(10))
    count = db.Column(db.Integer, default=0)
