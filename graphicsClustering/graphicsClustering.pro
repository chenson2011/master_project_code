#-------------------------------------------------
#
# Project created by QtCreator 2013-03-25T15:19:01
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = graphicsClustering
TEMPLATE = app


SOURCES += main.cpp\
        widget.cpp \
    node.cpp \
    edge.cpp \
    graph.cpp \
    kmean.cpp

HEADERS  += widget.h \
    node.h \
    edge.h \
    graph.h \
    kmean.h

FORMS    += widget.ui

INCLUDEPATH += D:\eigen

INCLUDEPATH += D:\boost_1_53_0
