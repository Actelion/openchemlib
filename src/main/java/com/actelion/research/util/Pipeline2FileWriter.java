package com.actelion.research.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Modest v. Korff
 * Idorsia Pharmaceuticals Ltd.
 * 17.06.2021 Start implementation
 **/
public class Pipeline2FileWriter<T> {


    private FileWriter fw;

    public Pipeline2FileWriter(Pipeline<T> pipe, File fiOut) throws IOException {
        fw = new FileWriter(fiOut);

        Runnable runWrite = new Runnable() {
            @Override
            public void run() {
                try {
                    boolean started=false;

                    while(!pipe.wereAllDataFetched()){
                        T t = pipe.pollData();
                        if(t == null) {
                            try {
                                Thread.sleep(250);
                            } catch (InterruptedException e) {
                                e.printStackTrace();
                            }
                            continue;
                        }

                        try {
                            if(started)
                                fw.write("\n");
                            fw.write(t.toString());
                            started=true;
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                } finally {
                    try {
                        fw.close();
                        System.out.println("Pipeline2FileWriter finished writing " + fiOut.getAbsolutePath() + ".");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
            }
        };

        ExecutorService executorServiceReadApprovedSymbols = Executors.newSingleThreadExecutor();
        executorServiceReadApprovedSymbols.submit(runWrite);
        executorServiceReadApprovedSymbols.shutdown();
    }



}
