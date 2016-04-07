set terminal pngcairo size 1024,768 enhanced font 'Verdana,10'
set output 'curves.png'

set nokey

set ylabel "cart speed, m/s"
set xlabel "time, s"
show xlabel
show ylabel

colors = "#55BF87 #22DE9E #E0C31A #C2096A #DB1BF5 #C5B743 #F75C6F #A4F780 #1532E9 #CC4E74 #339EAE #2A676C #611689 #7AEACD #399880 #DC47BE #6D9C4F #9B1A2A #494231 #07B903 #ADFE5A #DB867D #6187CA #104537 #74F2DF #58D562"

plot for [i=1:104:2] 'curves.dat' using 1:(column(i+1)) lc rgb word(colors, i%26) with lines

#w(x) = w0*(1-exp(alpha*x))
#w0 = 10
#alpha = -.1
#fit w(x) "curves.dat" using 1:2 via w0,alpha
#
#plot 'curves.dat' using 1:2  title "23.99V"  with lines lt rgb "#55BF87", \
#               '' using 1:3  title "23.05V"  with lines lt rgb "#22DE9E", \
#               '' using 1:4  title "22.11V"  with lines lt rgb "#E0C31A", \
#               '' using 1:5  title "21.17V"  with lines lt rgb "#C2096A", \
#               '' using 1:6  title "20.23V"  with lines lt rgb "#DB1BF5", \
#               '' using 1:7  title "19.29V"  with lines lt rgb "#C5B743", \
#               '' using 1:8  title "18.34V"  with lines lt rgb "#F75C6F", \
#               '' using 1:9  title "17.40V"  with lines lt rgb "#A4F780", \
#               '' using 1:10 title "16.46V"  with lines lt rgb "#1532E9", \
#               '' using 1:11 title "15.52V"  with lines lt rgb "#CC4E74", \
#               '' using 1:12 title "14.58V"  with lines lt rgb "#339EAE", \
#               '' using 1:13 title "13.64V"  with lines lt rgb "#2A676C", \
#               '' using 1:14 title "12.70V"  with lines lt rgb "#611689", \
#               '' using 1:15 title "11.76V"  with lines lt rgb "#7AEACD", \
#               '' using 1:16 title "10.82V"  with lines lt rgb "#399880", \
#               '' using 1:17 title "9.88V"   with lines lt rgb "#DC47BE", \
#               '' using 1:18 title "8.93V"   with lines lt rgb "#6D9C4F", \
#               '' using 1:19 title "7.99V"   with lines lt rgb "#9B1A2A", \
#               '' using 1:20 title "7.05V"   with lines lt rgb "#494231", \
#               '' using 1:21 title "6.11V"   with lines lt rgb "#07B903", \
#               '' using 1:22 title "5.17V"   with lines lt rgb "#ADFE5A", \
#               '' using 1:23 title "4.23V"   with lines lt rgb "#DB867D", \
#               '' using 1:24 title "3.29V"   with lines lt rgb "#6187CA", \
#               '' using 1:25 title "2.35V"   with lines lt rgb "#104537", \
#               '' using 1:26 title "1.41V"   with lines lt rgb "#74F2DF", \
#               '' using 1:27 title "0.47V"   with lines lt rgb "#58D562", \
# \
#               '' using 1:28 title "23.99V"  with lines lt rgb "#55BF87", \
#               '' using 1:29 title "23.05V"  with lines lt rgb "#22DE9E", \
#               '' using 1:30 title "22.11V"  with lines lt rgb "#E0C31A", \
#               '' using 1:31 title "21.17V"  with lines lt rgb "#C2096A", \
#               '' using 1:32 title "20.23V"  with lines lt rgb "#DB1BF5", \
#               '' using 1:33 title "19.29V"  with lines lt rgb "#C5B743", \
#               '' using 1:34 title "18.34V"  with lines lt rgb "#F75C6F", \
#               '' using 1:35 title "17.40V"  with lines lt rgb "#A4F780", \
#               '' using 1:36 title "16.46V"  with lines lt rgb "#1532E9", \
#               '' using 1:37 title "15.52V"  with lines lt rgb "#CC4E74", \
#               '' using 1:38 title "14.58V"  with lines lt rgb "#339EAE", \
#               '' using 1:39 title "13.64V"  with lines lt rgb "#2A676C", \
#               '' using 1:40 title "12.70V"  with lines lt rgb "#611689", \
#               '' using 1:41 title "11.76V"  with lines lt rgb "#7AEACD", \
#               '' using 1:42 title "10.82V"  with lines lt rgb "#399880", \
#               '' using 1:43 title "9.88V"   with lines lt rgb "#DC47BE", \
#               '' using 1:44 title "8.93V"   with lines lt rgb "#6D9C4F", \
#               '' using 1:45 title "7.99V"   with lines lt rgb "#9B1A2A", \
#               '' using 1:46 title "7.05V"   with lines lt rgb "#494231", \
#               '' using 1:47 title "6.11V"   with lines lt rgb "#07B903", \
#               '' using 1:48 title "5.17V"   with lines lt rgb "#ADFE5A", \
#               '' using 1:49 title "4.23V"   with lines lt rgb "#DB867D", \
#               '' using 1:50 title "3.29V"   with lines lt rgb "#6187CA", \
#               '' using 1:51 title "2.35V"   with lines lt rgb "#104537", \
#               '' using 1:52 title "1.41V"   with lines lt rgb "#74F2DF", \
#               '' using 1:53 title "0.47V"   with lines lt rgb "#58D562", \
# \
#               '' using 1:54 title "23.99V"  with lines lt rgb "#55BF87", \
#               '' using 1:55 title "23.05V"  with lines lt rgb "#22DE9E", \
#               '' using 1:56 title "22.11V"  with lines lt rgb "#E0C31A", \
#               '' using 1:57 title "21.17V"  with lines lt rgb "#C2096A", \
#               '' using 1:58 title "20.23V"  with lines lt rgb "#DB1BF5", \
#               '' using 1:59 title "19.29V"  with lines lt rgb "#C5B743", \
#               '' using 1:60 title "18.34V"  with lines lt rgb "#F75C6F", \
#               '' using 1:61 title "17.40V"  with lines lt rgb "#A4F780", \
#               '' using 1:62 title "16.46V"  with lines lt rgb "#1532E9", \
#               '' using 1:63 title "15.52V"  with lines lt rgb "#CC4E74", \
#               '' using 1:64 title "14.58V"  with lines lt rgb "#339EAE", \
#               '' using 1:65 title "13.64V"  with lines lt rgb "#2A676C", \
#               '' using 1:66 title "12.70V"  with lines lt rgb "#611689", \
#               '' using 1:67 title "11.76V"  with lines lt rgb "#7AEACD", \
#               '' using 1:68 title "10.82V"  with lines lt rgb "#399880", \
#               '' using 1:69 title "9.88V"   with lines lt rgb "#DC47BE", \
#               '' using 1:70 title "8.93V"   with lines lt rgb "#6D9C4F", \
#               '' using 1:71 title "7.99V"   with lines lt rgb "#9B1A2A", \
#               '' using 1:72 title "7.05V"   with lines lt rgb "#494231", \
#               '' using 1:73 title "6.11V"   with lines lt rgb "#07B903", \
#               '' using 1:74 title "5.17V"   with lines lt rgb "#ADFE5A", \
#               '' using 1:75 title "4.23V"   with lines lt rgb "#DB867D", \
#               '' using 1:76 title "3.29V"   with lines lt rgb "#6187CA", \
#               '' using 1:77 title "2.35V"   with lines lt rgb "#104537", \
#               '' using 1:78 title "1.41V"   with lines lt rgb "#74F2DF", \
#               '' using 1:79 title "0.47V"   with lines lt rgb "#58D562", \
# \
#               '' using 1:80 title "23.99V"  with lines lt rgb "#55BF87", \
#               '' using 1:81 title "23.05V"  with lines lt rgb "#22DE9E", \
#               '' using 1:82 title "22.11V"  with lines lt rgb "#E0C31A", \
#               '' using 1:83 title "21.17V"  with lines lt rgb "#C2096A", \
#               '' using 1:84 title "20.23V"  with lines lt rgb "#DB1BF5", \
#               '' using 1:85 title "19.29V"  with lines lt rgb "#C5B743", \
#               '' using 1:86 title "18.34V"  with lines lt rgb "#F75C6F", \
#               '' using 1:87 title "17.40V"  with lines lt rgb "#A4F780", \
#               '' using 1:88 title "16.46V"  with lines lt rgb "#1532E9", \
#               '' using 1:89 title "15.52V"  with lines lt rgb "#CC4E74", \
#               '' using 1:90 title "14.58V"  with lines lt rgb "#339EAE", \
#               '' using 1:91 title "13.64V"  with lines lt rgb "#2A676C", \
#               '' using 1:92 title "12.70V"  with lines lt rgb "#611689", \
#               '' using 1:93 title "11.76V"  with lines lt rgb "#7AEACD", \
#               '' using 1:94 title "10.82V"  with lines lt rgb "#399880", \
#               '' using 1:95 title "9.88V"   with lines lt rgb "#DC47BE", \
#               '' using 1:96 title "8.93V"   with lines lt rgb "#6D9C4F", \
#               '' using 1:97 title "7.99V"   with lines lt rgb "#9B1A2A", \
#               '' using 1:98 title "7.05V"   with lines lt rgb "#494231", \
#               '' using 1:99 title "6.11V"   with lines lt rgb "#07B903", \
#               '' using 1:100 title "5.17V"   with lines lt rgb "#ADFE5A", \
#               '' using 1:101 title "4.23V"   with lines lt rgb "#DB867D", \
#               '' using 1:102 title "3.29V"   with lines lt rgb "#6187CA", \
#               '' using 1:103 title "2.35V"   with lines lt rgb "#104537", \
#               '' using 1:104 title "1.41V"   with lines lt rgb "#74F2DF", \
#               '' using 1:105 title "0.47V"   with lines lt rgb "#58D562"
#
#

