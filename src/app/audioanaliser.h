#include <QObject>
#include <QQmlEngine>

namespace Minuet {
    class AudioAnaliser : public QObject {
        Q_OBJECT
        public:
            Q_INVOKABLE QString getNoteName(QString filename);
        signals:

        public slots:

    };
}
