SIM_FM_FILE=stripe_320x185_Co_fm-up_D26e-4_h25e-2nm.sh
SIM_SK_FILE=stripe_320x185_Co_sk-down_D26e-4_h25e-2nm.sh

# -----------------------------------------------------------------------------

image:
	docker build -t nebm .

image_clean:
	docker rm $(docker stop $(docker ps -q))
	docker rmi -f nebm

# -----------------------------------------------------------------------------

relaxation:
	cd ../sims/relaxation && \
		docker run -v `pwd`:/io -u $(id -u):$(id -g) \
		nebm /bin/bash -c "bash ${SIM_FM_FILE}"
	
	cd ../sims/relaxation && \
		docker run -v `pwd`:/io -u $(id -u):$(id -g) \
		nebm /bin/bash -c "bash ${SIM_SK_FILE}"
