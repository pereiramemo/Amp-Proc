SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-1.1   -t ghcr.io/epereira/amp-proc/module-1.1:latest   .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-1.2   -t ghcr.io/epereira/amp-proc/module-1.2:latest   .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-1.3   -t ghcr.io/epereira/amp-proc/module-1.3:latest   .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-2.1   -t ghcr.io/epereira/amp-proc/module-2.1:latest   .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-2.2.1 -t ghcr.io/epereira/amp-proc/module-2.2.1:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-2.2.2 -t ghcr.io/epereira/amp-proc/module-2.2.2:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-2.2.3 -t ghcr.io/epereira/amp-proc/module-2.2.3:latest .
docker build --network=host -f ${SCRIPT_DIR}/Dockerfile.module-3     -t ghcr.io/epereira/amp-proc/module-3:latest     .
