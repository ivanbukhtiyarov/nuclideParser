TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt


QMAKE_CXXFLAGS += -std=c++14
QMAKE_CXXFLAGS += -std=gnu++14
SOURCES += main.cpp \
    parse.cpp \
    nuclide.cpp \
    chain.cpp \
    pugixml.cpp \
    configure.cpp \
    functionals.cpp \
    materials.cpp \
    reactions.cpp \
    matrix.cpp \
    filter.cpp

HEADERS += \
    nuclide.h \
    parse.h \
    pugiconfig.h \
    chain.h \
    pugixml.h \
    pugiconfig.h \
    uncertainty.h \
    configure.h \
    functionals.h \
    materials.h \
    reactions.h \
    matrix.h \
    timeproc.h \
    filter.h
