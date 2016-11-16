# If it is not defined, we set the DMI value here
DMI?=32

SIM_FM_FILE=stripe_320x185_Co_fm-up_D${DMI}e-4_h25e-2nm.sh
SIM_SK_FILE=stripe_320x185_Co_sk-down_D${DMI}e-4_h25e-2nm.sh
NEBM_SIM_FOLDER=neb_stripe_320x185_Co_sk-down_fm-up_D${DMI}e-4_h25e-2nm_GEODESIC

# Climbing images
CI26=21
CI28=21
CI30=20
CI32=15
CI34=15
CI36=15

CI = $(CI${DMI})

# -----------------------------------------------------------------------------

# Build the image using the Dockerfile parameters
image:
	docker build -t nebm .

# Clean all the images
clean_images:
	@echo "This will remove the Docker nebm image "
	# docker rm `docker ps --filter=status=exited --filter=status=created -q`
	# docker rmi `docker images -a --filter=dangling=true -q`
	docker rmi -f nebm

# -----------------------------------------------------------------------------

# -t run a pseudo TTY, -d detach to run in the background
# -u is to run as user:group so files are not owned by root (somehow it is
# important to run the commands with -t) 
relaxation:
	@echo "Starting relaxation of states for D = ${DMI}e-4 Jm**-2"
	@echo "------------------------------------------------------"
	@cd ../sims/relaxation && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "bash ${SIM_FM_FILE}"
	@cd ../sims/relaxation && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "bash ${SIM_SK_FILE}"

nebm:
	@echo "Starting simulations for D = ${DMI}e-4 Jm**-2"
	@echo "---------------------------------------------"
	@cd ../sims && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "cd nebm/${NEBM_SIM_FOLDER} && bash generate_sk-disp_initial_state.sh"
	@cd ../sims && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "cd nebm/${NEBM_SIM_FOLDER} && bash initiate_simulation_skdisp_k1e4.sh"
	@cd ../sims && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "cd nebm/${NEBM_SIM_FOLDER} && bash initiate_simulation_k1e4.sh"
	@cd ../sims && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "cd nebm/${NEBM_SIM_FOLDER} && bash initiate_simulation_climbing-${CI}_k1e4.sh"

plot_nebm:
	@echo "Plotting data from NEBM simulation with D = ${DMI}e-4 Jm**-2"
	@echo "------------------------------------------------------------"
	@cd ../ && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "python plot/plot_energy_bands.py --method boundary --D_list ${DMI} --snapshots"
	@cd ../ && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "python plot/plot_energy_bands.py --method linear_interpolations --D_list ${DMI} --snapshots"
	@cd ../ && \
		docker run -t -v `pwd`:/io -u `id -u`:`id -g` \
		nebm /bin/bash -c "python plot/plot_energy_bands.py --method climbing --D_list ${DMI} --snapshots"
