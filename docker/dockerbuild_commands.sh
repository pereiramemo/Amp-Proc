SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

docker build --network=host -f ${SCRIPT_DIR}/1.1-quality-check.Dockerfile   -t ghcr.io/epereira/amp-proc/1.1-quality-check:latest   .
docker build --network=host -f ${SCRIPT_DIR}/1.2-primers-check.Dockerfile   -t ghcr.io/epereira/amp-proc/1.2-primers-check:latest   .
docker build --network=host -f ${SCRIPT_DIR}/1.3-primers-removal.Dockerfile   -t ghcr.io/epereira/amp-proc/1.3-primers-removal:latest   .
docker build --network=host -f ${SCRIPT_DIR}/2.1-dada2-pipeline.Dockerfile   -t ghcr.io/epereira/amp-proc/2.1-dada2-pipeline:latest   .
docker build --network=host -f ${SCRIPT_DIR}/2.2.1-vsearch-pipeline.Dockerfile -t ghcr.io/epereira/amp-proc/2.2.1-vsearch-pipeline:latest .
docker build --network=host -f ${SCRIPT_DIR}/2.2.2-vsearch-pipeline.Dockerfile -t ghcr.io/epereira/amp-proc/2.2.2-vsearch-pipeline:latest .
# docker build --network=host -f ${SCRIPT_DIR}/2.2.3-otu-to-seqtable.Dockerfile -t ghcr.io/epereira/amp-proc/2.2.3-otu-to-seqtable:latest .
# docker build --network=host -f ${SCRIPT_DIR}/3-taxa-annot.Dockerfile     -t ghcr.io/epereira/amp-proc/3-taxa-annot:latest     .