AM_CPPFLAGS = -I. -Iblas

noinst_LTLIBRARIES = liblinear.la
noinst_HEADERS = linear.h tron.h blas/blas.h

liblinear_la_SOURCES = \
linear.cpp \
tron.cpp \
blas/daxpy.c \
blas/ddot.c \
blas/dnrm2.c \
blas/dscal.c

noinst_PROGRAMS = train predict

train_SOURCES = train.cpp
train_LDADD = liblinear.la

predict_SOURCES = predict.cpp
predict_LDADD = liblinear.la


