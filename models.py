from flask_sqlalchemy import SQLAlchemy
from datetime import datetime

db = SQLAlchemy()

# --- MUNDO 1: APRENDIZADO (FLASHCARDS & LEITURA) ---
class Deck(db.Model):
    id = db.Column(db.String(50), primary_key=True)
    name = db.Column(db.String(100))
    icon = db.Column(db.String(10))
    cards = db.relationship('Card', backref='deck', lazy=True, cascade="all, delete-orphan")

class Card(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    front = db.Column(db.String(200)) 
    back = db.Column(db.String(200))
    ipa = db.Column(db.String(200)) # Fonética (IPA)
    context = db.Column(db.Text)
    deck_id = db.Column(db.String(50), db.ForeignKey('deck.id'))
    
    # SRS (Sistema de Repetição Espaçada)
    next_review = db.Column(db.String(20)) # YYYY-MM-DD
    interval = db.Column(db.Integer, default=0)
    ease_factor = db.Column(db.Float, default=2.5)

class Story(db.Model):
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

# --- MUNDO 2: PESQUISA ACADÊMICA (TESE & PROJETOS) ---
class Project(db.Model):
    id = db.Column(db.String(50), primary_key=True)
    title = db.Column(db.String(200))
    target_journal = db.Column(db.String(100))

class Reference(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    title = db.Column(db.String(300))
    authors = db.Column(db.String(200))
    year = db.Column(db.String(4))
    status = db.Column(db.String(20)) # 'to_read', 'done'
    pdf_filename = db.Column(db.String(200))
    abstract = db.Column(db.Text)
    project_id = db.Column(db.String(50), db.ForeignKey('project.id'))

class Draft(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    section_name = db.Column(db.String(50))
    content = db.Column(db.Text)
    last_updated = db.Column(db.DateTime, default=datetime.utcnow)
    project_id = db.Column(db.String(50), db.ForeignKey('project.id'))
