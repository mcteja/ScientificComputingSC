
if(trigger==1)
s = SplashScreen( 'Splashscreen', 'splashimage.png', ...
                        'ProgressBar', 'on', ...
                        'ProgressPosition', 5, ...
                        'ProgressRatio', 0.4 );
      %s.addText( 300, 300, 'My Cool App', 'FontSize', 30, 'Color', [0 0 0.6] )
%s.addText( 300, 240, 'v1.0', 'FontSize', 20, 'Color', [0.2 0.2 0.5] )
s.addText( 50, 300, 'Loading...', 'FontSize', 15, 'Color', 'white' )
end
if(trigger==2)
delete(s)
end
        
return