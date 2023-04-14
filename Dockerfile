FROM continuumio/miniconda3

LABEL author="https://github.com/a-hr"
LABEL description="Wrapper for UMIclusterer (https://github.com/a-hr/UMIclusterer). Includes all its dependencies."

WORKDIR /usr/bin/UMIclusterer

COPY ./ .

RUN conda install python=3.8.15

RUN conda install --file spec-file.txt

RUN pip install levenshtein

CMD [ "/bin/bash" ]