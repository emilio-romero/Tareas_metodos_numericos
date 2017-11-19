CFLAGS=-c -std=c99 -lm 
EXEC=ejecutable 

CPPS:=$(wildcard src/*.c)
OBJS:=$(patsubst src/%.c,obj/%.o,$(CPPS))
OBJDIR:=obj

all: $(EXEC)

$(EXEC): $(OBJS)
	gcc -o $@ $(OBJS) -std=c99 -lm

$(OBJDIR)/%.o: src/%.c | $(OBJDIR)
	gcc $(CFLAGS) $< -o $@ 
$(OBJDIR): 
	mkdir -p $@
clean: 
	rm -f -r $(OBJDIR) $(EXEC)
