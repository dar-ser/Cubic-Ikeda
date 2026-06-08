# ── project ──────────────────────────────────────────────────────────────────
PROGS   = generateImages # visual # periodicOrbit 
OTHERS  = computations/bifurcation computations/pseudospectral computations/newton computations/cover computations/attractor computations/pca

# ── CAPD ──────────────────────────────────────────────────────────────────────
CAPDBINDIR =
CAPDCXX   := $(shell $(CAPDBINDIR)capd-config --variable=capd_cxx)
CAPDFLAGS  = $(shell $(CAPDBINDIR)capd-config --cflags)
CAPDLIBS   = $(shell $(CAPDBINDIR)capd-config --libs)

# ── compiler ──────────────────────────────────────────────────────────────────
CXXFLAGS += $(CAPDFLAGS) -O2 -std=c++17 -Icomputations  -I/mnt/c/Users/darin/Documents/Studia/eigen-5.0.0

# ── directories ──────────────────────────────────────────────────────────────
OBJDIR = .obj/

# ── object files ─────────────────────────────────────────────────────────────
OTHERS_OBJ = $(OTHERS:%=$(OBJDIR)%.o)
OBJ_FILES  = $(OTHERS_OBJ) $(PROGS:%=$(OBJDIR)%.o)

.PHONY: all clean run

# ── build all programs ───────────────────────────────────────────────────────
all: $(PROGS:%=output/%)

# ── linking rule (CAPD standard) ─────────────────────────────────────────────
output/%: $(OBJDIR)%.o $(OTHERS_OBJ)
	@mkdir -p output
	$(CAPDCXX) -o $@ $< $(OTHERS_OBJ) $(CAPDLIBS)

# ── dependency tracking ──────────────────────────────────────────────────────
-include $(OBJ_FILES:%=%.d)

# ── compilation rule ──────────────────────────────────────────────────────────
$(OBJ_FILES): $(OBJDIR)%.o : %.cpp
	@mkdir -p $(dir $@)
	$(CAPDCXX) $(CXXFLAGS) -MT $@ -MD -MP -MF $(@:%=%.d) -c -o $@ $<

# ── run ──────────────────────────────────────────────────────────────────────
run: all
	./output/generateImages
# 	./output/periodicOrbit
# 	./output/visual

# ── clean ────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(OBJDIR) output