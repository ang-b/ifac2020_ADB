%% output
figure(1)

for i=1:N
    ax(i) = subplot(3,2,i);
    plot(t, reshape(yT(1,i,:), [1 nTot]), ...
         t, reshape(y(1,i,:), [1 nTot]), '--', ...
         t, xref(1,:) );
    %ylim([-.5 .5]);
    legend({sprintf('$\\tilde y_{%i}$', i), ...
            sprintf('$y_{%i}$', i), ...
            sprintf('$y_{%i, ref}$', i)} , ...
            'Interpreter', 'latex', 'Location', 'NorthWest');
    linkaxes(ax,'x');
end


% %% decoupled error
% figure(2)
% for i=1:N
%     subplot(3,2,i);
%     plot(t, reshape(yT(1,i,:) - xd(1,i,:), [1 nTot]));
% end
% 
% %% coupled error
% figure(3)
% for i=1:N
%     subplot(3,2,i);
%     plot(t, reshape(yT(1,i,:) - xc(1,i,:), [1 nTot]));
% end

%% residuals
figure(4)
for i=1:N
    subplot(3,2,i);
    plot(t, dcres(i,:).');
    ylabel(sprintf('$\\|r_{%i}\\| $', i), ...
            'Interpreter', 'latex', 'FontSize', 12);
    ylim([0, 0.05]);
end
% 
% %% attacker
% figure(5)
% plot(t, xA)
% 
% %% inputs
% figure(6)
% for i=1:N
%     subplot(3,2,i);
%     plot(t, reshape(u(:,i,:) - utilde(:,i,:), [mI nTot]));
% end
