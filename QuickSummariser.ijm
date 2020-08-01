path=getDirectory("");
files=getFileList(path);
first=true;
for (i=0; i<lengthOf(files); i++) {
	if (endsWith(files[i], "_trackstats-driftcorrected.txt")==true) {
		str=File.openAsString(path+files[i]);
		lines=split(str, "\r\n");
		data0=split(lines[0], " ");
		if (first==true) {
			str=""+"file"+" "+data0[0];
			for (k=1; k<lengthOf(data0); k++) {
				str+=" "+data0[k];
			}
			print(str);
			first=false;
		}
		average=newArray(lengthOf(data0));
		count=0;
		for (j=1; j<lengthOf(lines); j++) {
			data=split(lines[j], " ");
			if (lengthOf(data)>=lengthOf(data0)) {
				for (k=0; k<lengthOf(data0); k++) {
					data[k]=parseFloat(data[k]);
					average[k]=average[k]+data[k];
				}
				count++;
			}
		}
		str=""+files[i]+" ";
		for (k=0; k<lengthOf(data0); k++) {
			average[k]=average[k]/count;
			str+=" "+average[k];
		}
		print(str);
	}
}
