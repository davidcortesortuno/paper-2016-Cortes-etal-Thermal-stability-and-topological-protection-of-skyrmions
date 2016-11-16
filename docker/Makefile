SIM_FM_FILE=stripe_320x185_Co_fm-up_D26e-4_h25e-2nm.sh
SIM_SK_FILE=stripe_320x185_Co_sk-down_D26e-4_h25e-2nm.sh
NEBM_SIM_FOLDER=neb_stripe_320x185_Co_sk-down_fm-up_D26e-4_h25e-2nm_GEODESIC

# -----------------------------------------------------------------------------

# Build the image using the Dockerfile parameters
image:
	docker build -t nebm .

# Clean all the images
clean_images:
	echo "This will remove the Docker nebm image "
	# docker rm `docker ps --filter=status=exited --filter=status=created -q`
	# docker rmi `docker images -a --filter=dangling=true -q`
	docker rmi -f nebm

# -----------------------------------------------------------------------------

# -t run a pseudo TTY, -d detach to run in the background
# -u is to run as user:group so files are not owned by root (somehow it is
# important to run the commands with -t) 
relaxation:
	cd ../sims/relaxation && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "bash ${SIM_FM_FILE}"
	cd ../sims/relaxation && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "bash ${SIM_SK_FILE}"

nebm_linear_interpolations:
	cd ../sims && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "cd nebm/${NEBM_SIM_FOLDER} && bash initiate_simulation_k1e4.sh"
	cd ../sims && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "cd nebm/${NEBM_SIM_FOLDER} && bash initiate_simulation_climbing-23_k1e4.sh"

plot_linear_interpolations:
	cd ../ && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "python plot/plot_energy_bands.py"
