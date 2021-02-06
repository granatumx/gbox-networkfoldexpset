FROM granatumx/gbox-py-sdk:1.0.0

RUN pip install networkx

RUN apt-get install -y graphviz
RUN apt-get install -y libgraphviz-dev
RUN pip install pygraphviz

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
