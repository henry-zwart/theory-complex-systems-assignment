MAIN ?= main.typ
PORT ?= 23625

EXPORT_NAME ?= solutions.pdf
COMPILE_ARGS ?= --root "$$(pwd)" $(EXPORT_NAME)

NEURON_MODELLING_RESULTS = results/neuron_modelling.json results/figures/neuron_interspiking_distribution.pdf
REPORT_DEPENDENCIES = report/main.typ report/solutions.typ $(NEURON_MODELLING_RESULTS)

.PHONY: clean preview

all: $(EXPORT_NAME)

$(EXPORT_NAME): $(REPORT_DEPENDENCIES)
	typst compile report/main.typ $(COMPILE_ARGS)

preview:
	typst compile report/main.typ $(COMPILE_ARGS) && \
	zathura solutions.pdf & \
	typst watch report/main.typ $(COMPILE_ARGS)


# === Experimental results
# = Neuron modelling
results/neuron_modelling.json: scripts/model_neuron.py | results results/figures
	uv run $<

results/figures: | results
	mkdir $@

results:
	mkdir $@

clean:
	rm -rf $(EXPORT_NAME)
