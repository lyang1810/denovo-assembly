package assembly.data.loader;

import assembly.data.DenovoPeptide;

import java.io.File;
import java.util.List;

public interface ResultLoader {
    List<DenovoPeptide> load(File[] files);
}
