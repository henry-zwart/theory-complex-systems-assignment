MAIN ?= main.typ
PORT ?= 23625

EXPORT_NAME ?= solutions.pdf
COMPILE_ARGS ?= --root "$$(pwd)" $(EXPORT_NAME)

REPORT_DEPENDENCIES = report/main.typ report/solutions.typ

.PHONY: clean preview

all: $(EXPORT_NAME)

$(EXPORT_NAME): $(REPORT_DEPENDENCIES)
	typst compile report/main.typ $(COMPILE_ARGS)

preview:
	typst compile report/main.typ $(COMPILE_ARGS) && \
	zathura solutions.pdf & \
	typst watch report/main.typ $(COMPILE_ARGS)

clean:
	rm -rf $(EXPORT_NAME)
