%% output
figure(1)

for i=1:N
    ax(i) = subplot(3,2,i);
    grid on
    plot(t, reshape(yT(1,i,:), [1 nTot]), ...
         t, reshape(y(1,i,:), [1 nTot]), '--', ...
         t, xref(1,:) );
    %ylim([-.5 .5]);
    legend({sprintf('$\\tilde y_{%i}$', i), ...
            sprintf('$y_{%i}$', i), ...
            sprintf('$y_{%i, ref}$', i)} , ...
            'Interpreter', 'latex', 'Location', 'NorthWest');
    %linkaxes(ax,'x');
end


% %% decoupled error
% figure(2)
% for i=1:N
%     subplot(3,2,i);
%     plot(t, reshape(yT(1,i,:) - xd(1,i,:), [1 nTot]));
% end

%% residuals

figure(3)
for i=1:N
    ebi = max(abs( squeeze(yT(1,i,t > 5)) - xref(1, t > 5)' ));
    rbi = max(dcres(i,t > tA), [], 2);
    subplot(3,2,i);
    plot(t, dcres(i,:).', [t(1), t(end)], repmat(1.1*ebi, [1 2]), '--');
    ylabel(sprintf('$\\|r_{%i}\\| $', i), ...
            'Interpreter', 'latex', 'FontSize', 12);
    ylim([0, rbi*1.2]);
end
% 
% %% attacker
% figure(4)
% plot(t, xA)
% 

%% inputs
figure(5)
utildeol = reshape([utildec{nA,:}],[2, nTot-1]);
uol = reshape([ucell{nA,:}],[2, nTot-1]);

for i=1:N
    subplot(3,2,i);
    plot(t, reshape(u(:,i,:), [mI nTot]), ...
         t, reshape(utilde(:,i,:), [mI nTot]), '--');
    if i == nA
        hold on
        plot(t(1:end-1), utildeol.', '--', t(1:end-1), uol.')
        hold off
    end
end
