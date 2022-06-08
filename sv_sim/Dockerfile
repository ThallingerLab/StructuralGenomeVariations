FROM davidebolo1993/visor
# their are no tags for visor than latest :/

LABEL base.image="ubuntu:20.04"
MAINTAINER Veronika Schusterbauer

# RUN conda activate visorenv
RUN conda install -y -n visorenv numpy
RUN conda install -y -n visorenv pandas
# RUN conda install -y -n visorenv -c bioconda pyranges
RUN conda install -y -n visorenv -c conda-forge biopython
RUN conda install -y -n visorenv -c conda-forge pyvcf

RUN mkdir /scripts
RUN mkdir /data 

#RUN mkdir /art && cd /art && \
#	curl -L0 https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz > artbinmountrainier2016.06.05linux64.tgz && \
#	tar xzvf artbinmountrainier2016.06.05linux64.tgz && \
#	rm artbinmountrainier2016.06.05linux64.tgz
	
ENV DEBIAN_FRONTEND=noninteractive 	
RUN apt update && apt install -y art-nextgen-simulation-tools
RUN apt update && apt install -y ncbi-blast+

ADD bedmaker/bed_maker_6_1.py /scripts/
ADD bedmaker/bed_maker_supplementary.py /scripts/
ADD Data_analysis/main_20210301.py /scripts/
ADD scripts/* /scripts/
ADD data/* /data/


