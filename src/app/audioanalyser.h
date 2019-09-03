#include <QObject>
#include <QQmlEngine>

namespace Minuet {
    class AudioAnalyser : public QObject {
        Q_OBJECT
        public:
            explicit AudioAnalyser(QObject *parent = 0);
            QString analyzeFile(const QString& in);
            Q_INVOKABLE QString getNoteName(QString filename);
        signals:

        public slots:
            // QString analyzeFile(const QString& in);

    };
}
