FILES = rt3dv25.obj utlrtv25.obj adv30d.obj dsp30d.obj ssmrtv25.obj btnrtv25.obj gcg30d.obj rteqnv25.obj fmi30d.obj rtrctv25.obj solver123.obj solver45.obj rxns.lib jacrxns.lib
rt3dv25dll.exe : $(FILES)
        df /exe:rt3dv25dll.exe -fast $(FILES)
#Uncomment the following 4 lines if you need to rebuild the rxns.lib and jacrxns.lib files; only required if you do not have the *.lib files
#rxns.lib : rxns.f
#        df /LD rxns.f
#jacrxns.lib : jacrxns.f
#        df /LD jacrxns.f
rt3dv25.obj : rt3dv25.f
        df -fast -c rt3dv25.f
utlrtv25.obj : utlrtv25.f
        df -fast  -c utlrtv25.f
adv30d.obj : adv30d.f
        df -fast -c adv30d.f
dsp30d.obj : dsp30d.f
        df -fast -c dsp30d.f
ssmrtv25.obj : ssmrtv25.f
        df -fast -c ssmrtv25.f
gcg30d.obj : gcg30d.f
        df -fast -c gcg30d.f
btnrtv25.obj : btnrtv25.f
        df -fast -c btnrtv25.f
fmi30d.obj : fmi30d.f
        df -fast -c fmi30d.f
rtrctv25.obj : rtrctv25.f
        df -fast -c rtrctv25.f
rteqnv25.obj : rteqnv25.f
        df -fast -c rteqnv25.f
solver123.obj : solver123.f
        df -fast -c solver123.f
solver45.obj : solver45.f
        df -fast -c solver45.f
