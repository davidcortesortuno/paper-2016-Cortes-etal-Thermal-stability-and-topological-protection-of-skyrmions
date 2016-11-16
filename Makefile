example:
	@echo "Building Docker container:"
	@cd docker && make image
	@echo "Running NEBM simulations for D = 0.721 meV (D = 3.2 mJ m **-2)"
	@cd docker && make relaxation && make nebm && make plot_nebm

run_all:
	@echo "Building Docker container:"
	@cd docker && make image
	@echo "Running NEBM simulations for different DMI values"
	@cd docker && \
		for D in 26 28 30 32 34 36; do \
	   		export DMI=$$D && make relaxation && make nebm && make plot_nebm;  \
		done

clean_docker:
	@echo "Attempting to remove docker container: nebm"
	@cd docker && make clean_image
