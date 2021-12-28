require(ggplot2)

counts = read.table("read-counts.txt", sep="\t", header=T)
runtimes = read.table("runtimes.txt", sep="\t", header=T)

png(filename="runtimes.png", width=800, height=550, res=100)
ggplot(runtimes) +
  aes(x=tool, y=runtime_seconds, fill=tool) + 
  geom_col() + 
  scale_y_continuous(minor_breaks=seq(0, 1700, 100)) + 
  scale_fill_brewer(palette = "Paired") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.subtitle = element_text(hjust = 0.5)) + 
  labs(x="Tool", y="Runtime (seconds)", 
       title="Runtime of Counting Guides in 4 FASTQs from Sanson et al.",
       subtitle="(Smaller bars are better)")
dev.off()

png(filename="read-counts.png", width=800, height=550, res=100)
ggplot(counts) + 
  aes(fill=analysis, x=sample, y=matched_reads) +
  scale_fill_brewer(palette = "Paired") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.subtitle = element_text(hjust = 0.5)) + 
  geom_bar(position="dodge", stat="identity") + 
  labs(x="Sample", y="Reads Matched to Guides",
       title="Matched Reads in 4 FASTQs from Sanson et al.",
       subtitle="(Bigger bars are better)")  
dev.off()
