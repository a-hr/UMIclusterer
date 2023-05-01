FROM continuumio/miniconda3

LABEL author="https://github.com/a-hr"
LABEL description="Wrapper for UMIclusterer (https://github.com/a-hr/UMIclusterer). Includes all its dependencies."

WORKDIR /usr/local/bin/UMIclusterer

COPY ./ .

RUN conda install python=3.8.15

RUN conda install --file spec-file.txt

ENV PATH="/usr/local/bin/UMIclusterer:$PATH" 
RUN chmod +x /usr/local/bin/UMIclusterer/umiclusterer.py

CMD [ "/bin/bash" ]