function hhTransitions(N,dim)
    # Function to generate transition matrices for household model
    # Input: N is the household size

    # Initialize things
    Qinf = zeros(dim,dim);
    Qrec = zeros(dim,dim);
    Qext = zeros(dim,dim);
    Qwane = zeros(dim,dim);
    dataI = Array{Int64}(zeros(dim,3))
    m = 0;
    I = Array{Int64}(zeros(N+1,N+1))

    # To help remember where to store the variables
    for ss = 0:N
        for ii = 0:(N-ss)
            m = m + 1;
            I[ss+1,ii+1] = m
        end
    end

    # Describe the epidemiological transitions

    # Counter for susceptibles
    for ss = 0:N
        # Counter for infecteds
        for ii = 0:(N-ss)
            # If susceptibles and infecteds are more than 1, then infection within the household can occur
            if (ss > 0 && ii > 0)
                Qinf[I[ss+1,ii+1],I[ss,ii+2]] = ii*ss/(N-1);
            end

            # If infecteds are more than 1, recovery can occur
            if ii > 0
                # Rate of recovery
                Qrec[I[ss+1,ii+1],I[ss+1,ii]] = ii;
            end

            # For external infection - just keep track of susceptibles
            if ss > 0
                # Rate of within household infection
                Qext[I[ss+1,ii+1],I[ss,ii+2]] = ss;
            end

            # Loss of protection hence becoming susceptible again. Possible if N-ss-ii = rr > 0
            if (N-ss-ii) > 0
                # Rate of loss of protection
                Qwane[I[ss+1,ii+1],I[ss+2,ii+1]] = N-ss-ii;
            end

            # Store the relevant indices to help identify the household configurations
            dataI[I[ss+1,ii+1],:] = [ss, ii, N-ss-ii];
        end
    end

    Qinf = Qinf - diagm(vec(sum(Qinf,2)),0);
    Qrec = Qrec - diagm(vec(sum(Qrec,2)),0);
    Qext = Qext - diagm(vec(sum(Qext,2)),0);
    Qwane = Qwane - diagm(vec(sum(Qwane,2)),0);

    # Return
    return Qinf, Qrec, Qext, Qwane, dataI
end
