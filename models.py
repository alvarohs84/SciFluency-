from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

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
    next_review = db.Column(db.String(20))
    interval = db.Column(db.Integer, default=0)

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
