# ==============================================================================
# rash - Minimal Task Runner
# ==============================================================================

SHELL    := /bin/bash
GIT_ROOT := $(shell git rev-parse --show-toplevel 2>/dev/null)
ETC      := $(GIT_ROOT)/etc

CLS     := '\033[H\033[J'
cRESET  := '\033[0m'
cYELLOW := '\033[1;33m'

help: ## show help
	@awk 'BEGIN{FS=":.*##"} \
	      /^[a-zA-Z_%\/.~$$-]+:.*##/ \
	      {printf "  \033[36m%-20s\033[0m %s\n",$$1,$$2}' \
	      $(MAKEFILE_LIST)

pyclean: ## remove python temporaries
	@find $(GIT_ROOT) -type d \
	      \( -name __pycache__ -o -name .pytest_cache \
	         -o -name "*.egg-info" \) -exec rm -rf {} +

push: ## commit with prompted msg and push
	@read -p "Reason? " msg; \
	 git commit -am "$$msg"; git push; git status

etc/simplex5.csv: etc/simplex5_noisy.py
	@python3 $< > $@

.PHONY: eps_simplex
eps_simplex: etc/simplex5.csv ## sweep eps on noisy 5-simplex
	@$(MAKE) eps FILE=etc/simplex5.csv

lint: ## lint f=x.py
	@pylint --disable=C0103,C0114,C0116 $f

test: ## run all tests
	@python3 -B rash.py --all

sh: ## launch dev shell (banner + etc/bash.rc)
	@-echo -e $(CLS)$(cYELLOW); figlet -W -f slant rash; \
	  echo -e $(cRESET)
	@-bash --init-file $(ETC)/bash.rc -i

~/tmp/%.pdf: %.py $(MAKEFILE_LIST) ## .py ==> .pdf
	@mkdir -p ~/tmp
	@echo "pdf-ing $@ ... "
	@a2ps -Br --quiet --landscape --chars-per-line=65 \
	      --lines-per-page=100 --line-numbers=1 --borders=no \
	      --pro=color --columns=3 -M letter \
	      -o - $< | ps2pdf - $@
	@open $@

PYCCO := $(HOME)/tmp

$(PYCCO)/%.html : %.py
	@mkdir -p $(PYCCO)
	gawk -f $(GIT_ROOT)/etc/py2md.awk $< > $(PYCCO)/$<
	pycco -d $(PYCCO)  $(PYCCO)/$^
	echo "p {text-align: right;}" >> $(PYCCO)/pycco.css

eps: ## sweep eps; FILE=path optional
	@for e in 0.01 0.05 0.1 0.2 0.35 0.5; do \
		printf "eps=%s  " $$e; \
		python3 -B rash.py --eps=$$e $(if $(FILE),--file=$(FILE),) --dim; \
	done
