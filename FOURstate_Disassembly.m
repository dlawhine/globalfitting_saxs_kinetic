function M=FOURstate_Disassembly(t,param,matrix,c0,kfwd,kback,alphaS,alphaB,betaS,betaB)
    [tt,Mtemp] = ode15s(@(T,c)ode_FOUR_state_disass_model(T,c,kfwd,kback,...
    alphaS,alphaB,betaS,betaB),t,c0);
Mtemp = Mtemp';
M=Mtemp+matrix;
end