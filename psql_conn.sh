#!/bin/bash
docker run -it --rm --net host --name psql -v /tmp:/tmp postgres psql -h localhost -U postgres -d congelados
