# coding: utf-8
from sqlalchemy import Column, Float, ForeignKey, Integer, String, Text
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata


class Node(Base):
    __tablename__ = 'node'

    id = Column(Integer, primary_key=True)
    parent = Column(Integer)
    left = Column(Integer)
    right = Column(Integer)
    name = Column(String(20))
    length = Column(Float)
    height = Column(Float)


class Taxon(Base):
    __tablename__ = 'taxon'

    taxon_id = Column(Integer, primary_key=True)
    taxon = Column(Text)
    kingdom = Column(Text, nullable=False)
    _class = Column('class', Text, nullable=False, index=True)
    ord = Column(Text, nullable=False, index=True)
    family = Column(Text, nullable=False, index=True)
    genus = Column(Text, nullable=False, index=True)
    bin_uri = Column(Text, nullable=False, index=True)
    opentol_id = Column(Integer, index=True)


class Barcode(Base):
    __tablename__ = 'barcode'

    barcode_id = Column(Integer, primary_key=True)
    processid = Column(Integer, nullable=False, index=True)
    marker_code = Column(Text, nullable=False)
    nuc = Column(Text, nullable=False, index=True)
    country = Column(Text, index=True)
    taxon_id = Column(ForeignKey('taxon.taxon_id'), nullable=False, index=True)

    taxon = relationship('Taxon')
