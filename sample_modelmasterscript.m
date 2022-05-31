infilelocationpath = t.path.model + "\\Codes\\SIRage_macro_GP";

cd(infilelocationpath)

if re_simulate == 1
    dynare SIRage_macro
end


